#!/usr/bin/perl -w
use strict;
use Getopt::Std;

########## ########## Get paramter ########## ##########
my %opts;
getopts('hto:w:p:d:s:f:c:e:u:', \%opts);
$opts{c}=0.95 unless ($opts{c});#coincide SNP proportion when extending.
$opts{u}=0.55 unless ($opts{u});#phase 0 and 1 cutoff value of scoring
$opts{w} = 500 unless ($opts{w});#window size of seed region selection
$opts{p} = 0.25 unless ($opts{p});#upper heter snp cutoff, alt fre/seq depth
$opts{d} = 0.75 unless ($opts{d});#lowwer heter snp cutoff, alt fre/seq depth
$opts{e} = 0.3 unless ($opts{e});#cutoff value of seed (SNP) pattern selection

&help and &info and die if ($opts{h}); 
&help and die unless ($opts{o} && @ARGV);

system "mkdir -p $opts{o}/";
system "mkdir -p $opts{o}/TEMP";

open OUT , ">$opts{o}/TEMP/score.temp";
open HIT , ">$opts{o}/TEMP/hit.temp";#read extend marker hit
open MK , ">$opts{o}/TEMP/markerSNP.temp";#marker SNP
open PT , ">$opts{o}/TEMP/seedSelect.temp";
open PH , ">$opts{o}/TEMP/phaseSNP.temp";#phase 0 SNP
open LOG , ">$opts{o}/log.txt";#log file

print MK "Pos\tAlt\tFre\n";
print PH "Phase\tPos\tMaxAlt\n";
print PT "Pattern:Fre\n";
print OUT "Phase\tReadNO\tLen\tMarkerNum\tMarkerYes\tMarkerNo\tSeqError\tHomoSnp\n";
print HIT "Phase\tReadNO\tLen\tMarkerNum\tMarkerYes\tMarkerNo\tSeqError\tHomoSnp\n";

########## ########## Data structure ########## ##########
my $read= {#data format likly to language C struct format
    QNAME => [],
    START => [],
    END => [],
    LEN => [],
    SNPPOS => [],
    SNPALT => [],
    SNPQUAL => [],
};
my $step = 0;

########## ########## Detect all SNP in all read ########## ##########
my @bam = @ARGV;
&detectSnp("mismatch");
$step++;
print " - Finish ($step/16) BAM file SNP detection.\n";
print LOG " - Finish ($step/16) BAM file SNP detection.\n";

########## ########## Identify heterozygous SNP marker, seq error and homo SNP ########## ##########
my %filter;#filtered marker hash
my %seqError;#read snp that seem to be sequencing error
my %homoSnp;#read snp that seem to be homozygous snp

my $heterSnpNum = &traverseSnp(\%filter, \%seqError, \%homoSnp);
$step++;
print " - Finish ($step/16) heterozygous SNP markers: $heterSnpNum determination.\n";
print LOG " - Finish ($step/16) heterozygous SNP markers: $heterSnpNum determination.\n";

########## ########## Detect ref-allel SNP in all read ########## ##########
my %refSnp;#ref-allel snp {qname}{pos}
&detectSnp("match");
$step++;
print " - Finish ($step/16) ref-allel SNP detection.\n";
print LOG " - Finish ($step/16) ref-allel SNP detection.\n";

########## ########## Set seed ########## ##########
my $bestWin = &traverseForSeedRegion(\%filter);
$step++;
print " - Finish ($step/16) seed windows region setting. Window: $bestWin.\n";
print LOG " - Finish ($step/16) seed windows region setting. Window: $bestWin.\n";

#seed arrays setting
my @seed0;
my @seed1;
&seedSelect($bestWin, \%filter, \%seqError, \%homoSnp, \@seed0, \@seed1);#return 2 seed array in diff hap as seed.
print "-- Phase 0 seeds array (NO): @seed0.\n";
print "-- Phase 1 seeds array (NO): @seed1.\n";
print LOG "-- Phase 0 seeds array (NO): @seed0.\n";
print LOG "-- Phase 1 seeds array (NO): @seed1.\n";
$step++;
print " - Finish ($step/16) candidate seeds array selection.\n";
print LOG " - Finish ($step/16) candidate seeds array selection.\n";

#artificial seed or say region seed setting
my %seed0Snp;
my %seed1Snp;
my $artMark = &artSeed(\@seed0, \%seed0Snp, \%filter);
&artSeed(\@seed1, \%seed1Snp, \%filter);
print " - Finish ($step/16) two artificial region seeds setting; each seed contain heter-SNP marker: $artMark.\n";
print LOG " - Finish ($step/16) two artificial region seeds setting.\n";

#phasing
&Phase(\%seed0Snp, "phase.0", 0);
&Phase(\%seed1Snp, "phase.1", 1);

sub Phase{
    my ($X, $fileName, $pha) = @_;
    my %seedSnp = %$X;
    my $allRead = $#{$read->{LEN}} + 1;
    ########## ########## Initialize phase_0 SNP ########## ##########
    my %phaseSnp;

    &phaseInitial(\%phaseSnp, \%filter, \%refSnp);#initial heter snp marker pos to phase0
    &seedInitial(\%phaseSnp, \%filter, \%refSnp, $X);#initial seed to phase 0
    $step++;
    print " - Finish ($step/16) phase $pha SNPs initialization.\n";
    print LOG " - Finish ($step/16) phase $pha SNPs initialization.\n";

    ########## ########## Grow SNP tree (phase 0 SNP markers) ########## ##########
    my @range;#detect range of genome [0]:start pos, [1]:end pos
    $range[0] = ($bestWin-1)*$opts{w};
    $range[1] = ($bestWin+1)*$opts{w};
    my ($extendTime, $extendLen) = &readExtend(\@range, \%phaseSnp, \%filter, \%refSnp, \%seqError, \%homoSnp, $pha, $opts{c});
    $step++;
    print "-- Phase $pha seed extend length (bp): $extendLen extend times: $extendTime (all read: $allRead)\n";
    print LOG "-- Phase $pha seed extend length (bp): $extendLen extend times: $extendTime (all read: $allRead)\n";
    print PT "-- Phase $pha seed extend length (bp): $extendLen extend times: $extendTime (all read: $allRead)\n";
    print " - Finish ($step/16) phase $pha heter-SNP-marker tree growth.\n";

    ########## ########## Filter phase_0 heter SNP markers ########## ##########
    my %phaseSnpFilter;
    &markerFilter(\%phaseSnp, \%phaseSnpFilter, $pha);
    $step++;
    print " - Finish ($step/16) phase $pha heter SNP markers filtering.\n";
    print LOG " - Finish ($step/16) phase $pha heter SNP markers filtering.\n";

    ########## ########## Scoring all reads ########## ##########
    my %qnameMark;#hash of qname and markYes(hit marker) value
    my @markerVal;#arrary of marker Yes(hit) value for midium caculation
    &scoring(\%phaseSnpFilter, \@markerVal, \%qnameMark, \%filter, \%seqError, \%homoSnp, \%refSnp, $pha);
    $step++;
    print " - Finish ($step/16) all reads scoring.\n";
    print LOG " - Finish ($step/16) all reads scoring.\n";

    ########## ########## Determine 2 haplotype ########## ##########
    &printResult(\@markerVal, \%qnameMark, "$fileName.qname");
    $step++;
    print " - Finish ($step/16) phase $pha result printing.\n";
    print LOG " - Finish ($step/16) phase $pha result printing.\n";
}

close OUT;
close PH;
close MK;
close HIT;
system"rm -r $opts{o}/TEMP/" if ($opts{t});

########## ########## Veen check ########## ##########
&veen;


########## ########## Functions ########## ##########
sub detectSnp {
    my $l=0;#line number or say read number
    for my $f (sort @bam){
        #best hit read only
        my %bestHit;#to keep line number
        open IN , "samtools view $f |" ;
        while (<IN>){
            next if (/^@/);
            $l++;
            my ($qname)=(split /\t/)[0];
            $_=~/AS:i:(.*?)\s/;
            my $as=$1;#AS field
            $bestHit{$qname}{$as}="$_";
        }
        close IN;
        
        for my $kq (sort keys %bestHit){#kq qname key
            my @AS=reverse sort {$a<=>$b} keys %{$bestHit{$kq}};
        
            my %seq_hash;#sequence based (start from 1) hash {ref_pos}{base N.O}{ref}{alt}
            my @line_snp_pos;
            my @line_snp_alt;
            my @line_snp_qual;#qualty of each snp
        
            my ($qname,$chr,$pos,$cigar,$seq,$qual)=(split /\t/,$bestHit{$kq}{$AS[0]})[0,2,3,5,9,10];
            $bestHit{$kq}{$AS[0]}=~/MD:Z:(.*?)\t/;
            my $md="$1";
            print LOG "No MD field was found in SAM/BAM file." and die unless($md);
            my $now_ci=$cigar;
            my $now_md=$md;
            my $now_pos=$pos;#position in reference genome
            my $now_seq=$seq;
            my $now_base=0;#N.O. of base in sequence
            my $start = $pos;
            my $len = length($seq);#lenthgh of read
            my $end = $pos+length($seq);#hard clip base is not avalible
        
            #use cigar field to figure out base N.O, base, type.
            $now_pos=$pos-1;
            while($now_ci=~ /^([0-9]+)([A-Z])(.*)/){
                my $num=$1;
                my $type=$2;
                $now_ci=$3;
                if ($type eq "M"){
                    for my $i (1 .. $num){
                        $now_base++;
                        $now_pos++;
                        my $base=substr($now_seq,0,1,"");#base of sequenced seq
                        my $bQua=substr($qual,0,1,"");#base quality
                        $seq_hash{$now_pos}{$now_base}{$base}{$bQua}{"M"}="mis/match";
                    }
                }elsif ($type eq "I"){
                    for my $i (1 .. $num){
                        $now_base++;
                        $now_pos+=0;
                        my $base=substr($now_seq,0,1,"");
                        my $bQua=substr($qual,0,1,"");#base quality
                        $seq_hash{$now_pos}{$now_base}{$base}{$bQua}{"I"}="insert";
                    }
                }elsif ($type eq "D"){
                    for my $i (1 .. $num){
                        $now_base+=0;
                        $now_pos++;
                        my $bQua = "NA";
                        $seq_hash{$now_pos}{$now_base}{"NA"}{$bQua}{"D"}="delet";
                    }
                }elsif ($type eq "S"){
                    for my $i (1 .. $num){
                        $now_base++;
                        $now_pos+=0;
                        my $base=substr($now_seq,0,1,"");#base of sequenced seq
                        my $bQua=substr($qual,0,1,"");#base quality
                        $seq_hash{$now_pos}{$now_base}{"NA"}{$bQua}{"S"}="soft";
                    }
                }elsif ($type eq "H"){
                }else{
                    print LOG "Unknow Cigar character.\n";
                    die;
                }
            }
            #use MD field to figure out reference sequence
            my $krp=$pos-1;
            while( $now_md ){
                if ($now_md =~ /^([0-9]+)(.*)/){
                    my $cat=$1;
                    $now_md=$2;
                    if ($cat == 0){
                        next;
                    }else{
                        for my $i (1 .. $cat){
                            $krp++;
                            for my $ksp (keys %{$seq_hash{$krp}}){
                                for my $kb (keys %{$seq_hash{$krp}{$ksp}}){
                                    for my $kqu (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                                        for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}{$kqu}}){
                                            $seq_hash{$krp}{$ksp}{$kb}{$kqu}{$kt}="$kb" if ($kt eq "M");
                                            $seq_hash{$krp}{$ksp}{$kb}{$kqu}{$kt}="NA" if ($kt eq "I");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if ($now_md =~ /^([A-Z]+)(.*)/){
                    #insertion or mis-match
                    my $cat=$1;
                    $now_md=$2;
                    for my $i ( 1 .. length($cat) ){
                        my $ref=substr($cat,0,1,"");
                        $krp++;
                        for my $ksp (keys %{$seq_hash{$krp}}){
                            for my $kb (keys %{$seq_hash{$krp}{$ksp}}){
                                for my $kqu (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                                    for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}{$kqu}}){
                                        if ($kt eq "M"){
                                            delete $seq_hash{$krp}{$ksp}{$kb}{$kqu}{"M"};
                                            $seq_hash{$krp}{$ksp}{$kb}{$kqu}{"Mis"}="$ref"; 
                                        }
                                        print LOG "Script error:$ksp $krp $kt == M\n" and die if($kt eq "I");
                                    }
                                }
                            }
                        }
                    }
                }
                if ($now_md =~ /^\^([A-Z]+)(.*)/){
                    my $cat=$1;
                    $now_md=$2;
                    for my $i ( 1 .. length($cat) ){
                        $krp++;
                    }
                }
            }
            #check result
            for my $krp (sort {$a<=>$b} keys %seq_hash){
                for my $ksp (keys %{$seq_hash{$krp}}){
                    for my $kb (keys %{$seq_hash{$krp}{$ksp}}){
                        for my $kqu (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                            for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}{$kqu}}){
                                my $kr=$seq_hash{$krp}{$ksp}{$kb}{$kqu}{$kt};
                                print LOG "Script error:$krp\t$ksp\t$kb\t$kt\t$kr\n" and die if ($kt eq "Mis" && $kb eq $kr);
                                ########## important ########## 
                                if ($kt eq "Mis" && $_[0] eq "mismatch"){
                                    push @line_snp_pos , $krp; 
                                    push @line_snp_alt , $kb; 
                                    push @line_snp_qual , (ord($kqu)-33); 
                                }
                                elsif ($kt eq "M"  && $_[0] eq "match"){
                                    $refSnp{$qname}{$krp}=(ord($kqu)-33) if (exists $filter{$krp});
                                }
                            }
                        }
                    }
                }
            }
            if ($_[0] eq "mismatch"){
                push @{ $read->{QNAME} } , $qname;
                push @{ $read->{START} } , $start;
                push @{ $read->{END} } , $end;
                push @{ $read->{LEN} } , $len;
                push @{ $read->{SNPPOS} } , [ @line_snp_pos ];
                push @{ $read->{SNPALT} } , [ @line_snp_alt ];
                push @{ $read->{SNPQUAL} } , [ @line_snp_qual ];
            }
        }
    }
}

sub range_overlap{
    my $start=$_[0];#range start pos
    my $end=$_[1];#range end pos
    my $tarS=$_[2];#target read start pos
    my $tarE=$_[3];#target read end pos
    #print"$start $end $tarS $tarE\n";
    if ($tarE<$start){#read is out of range(left)
        return 0;
    }elsif($tarS>$end){#read is out of range(rigth)
        return 0;
    }else{#read has overlap with range
        return (($tarE<$end?$tarE:$end)-($tarS>$start?$tarS:$start));
    }
}

sub range_marker{
    my ($read_s , $read_e) = @_;#read start pos and end pos
    my $s = 0;#counter
    for my $kpos (keys %filter){
        if($kpos>=$read_s && $kpos<=$read_e){
            $s++;
        }
    }
    #print"$read_s $read_e ",$read_e-$read_s," $s\n";
    return $s;
}

sub max_alt{
    my $phaseSnp = shift;
    my $kpos = shift;

    my $maxAlt = "NA";
    my $max = 0;
    for my $kalt (keys %{$$phaseSnp{$kpos}}){
        if ($$phaseSnp{$kpos}{$kalt} > $max){
        $max = $$phaseSnp{$kpos}{$kalt};
        $maxAlt = $kalt;
        }
    }
    return $maxAlt;
}

sub traverseSnp{#traverse all reads to figure out heter-SNP
    my ($filter, $seqError, $homoSnp) = @_;
    my %depth;#seq-depth
    my %sAlt;#alt and frequence of every SNP
    my $heterSnpNum = 0;
    for my $i (0 .. $#{$read->{QNAME}}){
        for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            my $qual = $read->{SNPQUAL}[$i][$j];
            $sAlt{$pos}{$alt}++;
        }
        for my $z ($read->{START}[$i] .. $read->{END}[$i]){
            $depth{$z}++;
        }
    }
    for my $kpos (sort {$a<=>$b} keys %sAlt){
        for my $kalt (sort keys %{$sAlt{$kpos}}){
            my $rate = $sAlt{$kpos}{$kalt}/$depth{$kpos};
            if( $rate > $opts{p} && $rate < $opts{d}){
                $$filter{$kpos}{$kalt}++;
                $heterSnpNum++;
                print MK "$kpos\t$kalt\t$sAlt{$kpos}{$kalt}\n";
            }elsif($rate <= $opts{p}){
                $$seqError{$kpos}{$kalt}++;
            }elsif($rate >= $opts{d}){
                $$homoSnp{$kpos}{$kalt}++;
            }
        }
    }
    return $heterSnpNum;
}

sub seedInitial{
    my ($X, $Y, $Z, $W) = @_;
    my $phaseSnp = $X;
    my %filter = %$Y; 
    my %refSnp = %$Z;
    my %seedSnp = %$W;
    for my $pos (keys %seedSnp){
        my $alt = $seedSnp{$pos};
        $$phaseSnp{$pos}{$alt} += 1;
    }
}

sub phaseInitial{
    my ($x,$y,$z)= @_;
    my $phaseSnp = $x;
    my %filter = %$y;
    my %refSnp = %$z;
    for my $i (0 .. $#{$read->{QNAME}}){#initial phase 0 SNP of seed read
        for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            my $qual = $read->{SNPQUAL}[$i][$j];
            if (exists $filter{$pos}){#same to filtered marker
                $$phaseSnp{$pos}{$alt} += 0;
            }
        }
        for my $kpos ($read->{START}[$i] .. $read->{END}[$i]){
            if (exists $filter{$kpos}){
                if (exists $refSnp{$read->{QNAME}[$i]}{$kpos}){
                    $$phaseSnp{$kpos}{"R"} += 0;
                }
            }
        }
    }
}

sub readExtend{
    my ($A, $B, $C, $D, $E, $F, $pha, $cor) = @_;
    my @range = @$A;
    my $phaseSnp = $B;
    my %filter = %$C;
    my %refSnp = %$D;
    my %seqError = %$E;
    my %homoSnp = %$F;

    my $p0r = 0;#read number determined as phase_0
    my $pre_p0r = -1;#previous p0r
    my %readConsider;# read been considered
    my $ovl = 1; 
    my $rC = 1;#read been considered number
    for (my $maxRedo=1+$#{$read->{LEN}}; $maxRedo; $maxRedo--){
        my %accordantRead;#reads can be add to "detect region" {consistency}{i}
        last if ($ovl < 0.1);#the overlap is too low
        $rC = keys %readConsider;
        last if ($rC == (1+$#{$read->{LEN}}));
        print "range(bp): ",$range[1]-$range[0],"\tread consider: $rC\tread extend: $p0r","\toverlap: ",$ovl,"\n";
        for my $i (0 .. $#{$read->{QNAME}}){
            if (exists $readConsider{$i}) { # skip phased read
                next;
            }
            else{
                next if (&range_overlap($range[0],$range[1],$read->{START}[$i],$read->{END}[$i])<=($ovl*$read->{LEN}[$i]));
                    $readConsider{$i}++;#
                my %conReadSnp0;
                for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
                    my $pos = $read->{SNPPOS}[$i][$j];
                    my $alt = $read->{SNPALT}[$i][$j];
                    my $qual = $read->{SNPQUAL}[$i][$j];
                    $conReadSnp0{$i}{$pos}=$alt;
                }
                my $xy = &markValue($i, \%conReadSnp0, \%refSnp, \%filter, \@range, $B);#%phaseSnp=%$B
                if($xy >= $cor){#phase0
                    $accordantRead{$xy}{$i}++;#reads can be add to "detect region" {consistency}{i}
                }else{
                    #$readConsider{$i}++;#
                    #$rC++;
                }
            }
        }
        for my $kxy (reverse sort {$a<=>$b} keys %accordantRead){#reads can be add to "detect region"
            for my $ki (sort keys %{$accordantRead{$kxy}}){
                $p0r += &extendBestRead($ki, \%refSnp, \@range, $B);#@range = @$A; $phaseSnp = $B;
                $readConsider{$ki}++;#
                $rC++;
            }
        }
        if ($pre_p0r == $p0r){
            $ovl -= 0.01;
        }
        $pre_p0r = $p0r;
    }
    return $p0r, ($range[1]-$range[0]);
}

sub markValue{#belong to readExtend part
    my ($i, $X, $Y, $Z, $R, $phaseSnp) = @_;
    my %readSnp0 = %$X;
    my %refSnp = %$Y;
    my %filter = %$Z;
    my @range = @$R;
    my $numMark = &range_marker(($read->{START}[$i]>$range[0]?$read->{START}[$i]:$range[0]), ($read->{END}[$i]<$range[1]?$read->{END}[$i]:$range[1]));# marker alt should be detect
    my $markYes = 0;#SNP is same to phase_0 marker alt (position and alt-allel)
    my $markNo = 0;#SNP is different to phase_0 marker alt (ref-allel or other base)
    my $sE = 0; #sequencing error
    my $hS = 0; #homo snp 
    my $otherSnp = 0;#SNP is not exist in marker alt
    #my $rS = 0;#SNP is heterzogous one is snp one is ref
    for my $ki (keys %readSnp0){
        for my $pos (keys %{$readSnp0{$ki}}){
            my $alt = $readSnp0{$ki}{$pos};
            if ($pos>=$range[0] && $pos <=$range[1]){#only consider the overlap part SNP
                if(exists $filter{$pos}){
                    if ($alt eq &max_alt($phaseSnp, $pos)){#pos and alt is same(alt-allel)
                        $markYes++;
                    }
                    else{#ref-allel
                        $markNo++;
                    }
                }
                elsif(exists $seqError{$pos}{$alt}){
                    $sE++;
                }
                elsif(exists $homoSnp{$pos}{$alt}){
                    $hS++;
                }
                else{
                    $otherSnp++;
                }
            }
        }
    }
    for my $pos (($read->{START}[$i]>$range[0]?$read->{START}[$i]:$range[0]) .. ($read->{END}[$i]<$range[1]?$read->{END}[$i]:$range[1])){
        if (exists $refSnp{$read->{QNAME}[$i]}{$pos}){
            my $alt = "R";
            if ($alt eq &max_alt($phaseSnp, $pos)){#pos and alt is same(alt-allel)
                $markYes++;
            }
            else{#ref-allel
                $markNo++;
            }
        }
    }
    my $xy = 0;#($markYes/($markYes+$markNo))
    if($markYes+$markNo){
        $xy = $markYes/($markYes+$markNo);
    }
    return $xy;
}


sub extendBestRead{#after traverse all read, add the highest read only
    my ($X, $Y, $Z, $W) = @_;
    my $ki = $X;
    my $refSnp = %$Y;
    my $range = $Z;
    my $phaseSnp = $W;

    $$range[0]=($$range[0]<$read->{START}[$ki]?$$range[0]:$read->{START}[$ki]);#re-new start of range
    $$range[1]=($$range[1]>$read->{END}[$ki]?$$range[1]:$read->{END}[$ki]);#re-new end of range
    #print HIT "$pha\t$ki\t$read->{LEN}[$ki]\t$numMark\t$markYes\t$markNo\t$sE\t$hS\n";
    my %readSnp;
    my %readSnp0;
    for my $j (0 ..$#{$read->{SNPPOS}[$ki]}){
        my $pos = $read->{SNPPOS}[$ki][$j];
        my $alt = $read->{SNPALT}[$ki][$j];
        my $qual = $read->{SNPQUAL}[$ki][$j];
        $readSnp{$pos}{$alt}=$qual;
        $readSnp0{$pos}=$alt;
    }
    for my $z ($read->{START}[$ki] .. $read->{END}[$ki]){#z is position from read start to end
        if (exists $filter{$z}){#there should be SNP marker
            if (exists $readSnp{$z}){#there is read SNP
                for my $kalt (keys %{$readSnp{$z}}){
                    #SNP base qual -10log10{Pr}
                    my $kqual = 1-10**($readSnp{$z}{$kalt}/(0-10));#base qual of SNP (<1)
                    $$phaseSnp{$z}{$kalt}+=$kqual;
                }
            }
            elsif (exists $refSnp{$read->{QNAME}[$ki]}{$z}){#this read is ref-allel
                $$phaseSnp{$z}{"R"}+=1-10**(($refSnp{$read->{QNAME}[$ki]}{$z})/(0-10));
            }
            else{#complex situation
            #indels
            }
        }
    }
    return 1;
}

sub markerFilter{
    my ($x, $y, $pha) = @_;
    my %phaseSnp = %$x;
    my $phaseSnpFilter = $y;
    for my $kpos (keys %phaseSnp){
        my $max = 0.1;
        my $maxAlt = "NA";
        for my $kalt (sort keys %{$phaseSnp{$kpos}}){
            if ($phaseSnp{$kpos}{$kalt} > $max){
                $max = $phaseSnp{$kpos}{$kalt};
                $maxAlt = $kalt;
            }
        }
        $$phaseSnpFilter{$kpos}="$maxAlt";
        print PH "$pha\t$kpos\t$maxAlt\n";
    }
}

sub scoring{
    my ($A, $B, $C, $D, $E, $F, $G, $pha) = @_;
    my %phaseSnpFilter = %$A;
    my $markerVal = $B;
    my $qnameMark = $C;
    my %filter = %$D;
    my %seqError = %$E;
    my %homoSnp = %$F;
    my %refSnp = %$G;

    for my $i (0 .. $#{$read->{QNAME}}){
        my $numMark = &range_marker($read->{START}[$i], $read->{END}[$i]);
        my $markYes = 0;
        my $markNo = 0;
        my $otherSnp = 0;
        my $sE = 0;
        my $hS = 0;
        my %readSnp;
        my %readSnp0;
        for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            my $qual = $read->{SNPQUAL}[$i][$j];
            $readSnp{$pos}{$alt}=$qual;
            $readSnp0{$pos}="$alt";
        }
        for my $z ($read->{START}[$i] .. $read->{END}[$i]){
            if (exists $filter{$z}){#there should be SNP marker
                if (exists $readSnp{$z}){#there is read SNP
                    my $alt = $readSnp0{$z};
                    my $kqual = 1-10**($readSnp{$z}{$alt}/(0-10));#base qual of SNP (<1)
                    if ($alt eq $phaseSnpFilter{$z}){#read snp alt is same to phase0 library
                        $markYes+=$kqual;
                    }
                    elsif(exists $seqError{$z}{$alt}){#read snp is sequencing error
                        $sE++;
                    }
                    elsif(exists $homoSnp{$z}{$alt}){#read snp is homo-snp which is not marker
                        $hS++;
                    }
                }
                elsif($phaseSnpFilter{$z} eq "R"){#there is no read SNP because SNP marker is ref-alt
                    if (exists $refSnp{$read->{QNAME}[$i]}{$z}){#this read is ref-allel
                        $markYes+=1-10**(($refSnp{$read->{QNAME}[$i]}{$z})/(0-10));
                    }else{
                        $markNo++;
                    }
                }
                else{#differ to phase0 library
                    $markNo++;
                }
            }
        }
        print OUT "$pha\t$i\t$read->{LEN}[$i]\t$numMark\t$markYes\t$markNo\t$sE\t$hS\n";
        if ($numMark){
            $markYes=$markYes/$numMark; 
            $markNo=$markNo/$numMark;
        }else{
            #$markYes=1;#if there is no heter SNP, which means can not judge, consider this read belongs to all 2 hap 
            $markYes=0;#if there is no heter SNP, which means can not judge, consider this read belongs to no hap 
        }
        push @$markerVal , $markYes;
        $$qnameMark{$read->{QNAME}[$i]}=$markYes;
    }
}

sub printResult{
    my ($A, $B, $fileName) = @_;
    my @markerVal = @$A;
    my %qnameMark = %$B;

    my $count = @markerVal;
    my $p0=0;#reads number of phase_0
    my $p1=0;#reads number of phase_1
    my @sortMarkerVal =  sort @markerVal;#from small to large
    my $markerCutOff0 = $sortMarkerVal[int($opts{u}*($count))];#small
    print "-- Marker hit cutoff value: $markerCutOff0\n";
    #####;
    open OUT0 , ">$opts{o}/$fileName";
    for my $kq (keys %qnameMark){
        my $markYes = $qnameMark{$kq};
        if ($markYes > $markerCutOff0){#phase0
            print OUT0"$kq\t$fileName\n";
            $p0++;
        }
    }
    close OUT0;
}

sub traverseForSeedRegion{#traverse all reads to figure out high heter and higt coverage region
    my ($x) = @_;
    my %filter = %$x; 
    my %seqDepth;#sequencing depth of each window
    my %heterNum;#heter-snp-marker num of each window
    for my $i (0 .. $#{$read->{QNAME}}){
        my $start = $read->{START}[$i];
        my $end = $read->{END}[$i];
        my $s = int($start/$opts{w})+1;
        my $sr = ($start % $opts{w})/$opts{w};
        my $e = int($end/$opts{w});
        my $er = ($end % $opts{w})/$opts{w};
        for my $w ($s .. $e){
            $seqDepth{$w}++;
        }
        my $w = $s - 1;
        $seqDepth{$w} += $sr;
        $w = $e + 1;
        $seqDepth{$w} += $er;
    }
    for my $kpos (keys %filter){
        my $w = int($kpos/$opts{w})+1;
        $heterNum{$w}++;
    }
    #judice
    #high seq depth and high heterozygosity region selection
    my $maxSeqDepth = 0;
    for my $k (keys %seqDepth){
        $maxSeqDepth = $seqDepth{$k} if ($seqDepth{$k} > $maxSeqDepth);
    }
    my $maxHeterNum = 0;
    my $maxHNWin;
    for my $kwin (keys %heterNum){
        next if ($seqDepth{$kwin} < 0.4*$maxSeqDepth);
        next if (!exists $heterNum{($kwin+1)} || !exists $heterNum{($kwin-1)});
        my $tri = $heterNum{$kwin}+$heterNum{($kwin+1)}+$heterNum{($kwin-1)};#three neighbor windows heter SNP num
        if ($tri > $maxHeterNum){
            $maxHeterNum = $tri;
            $maxHNWin = $kwin;
            #print "heterNum\t$kwin\t$tri\n";
        }
    }
    return $maxHNWin;
}

sub seedSelect{
    my ($bestWin, $x, $y, $z, $s0, $s1) = @_;#$s0 is \@seed0 $s1 is \@seed1;
    my %filter = %$x; 
    my %seqError = %$y; 
    my %homoSnp = %$z;
    #find reads in the "region"
    my @seedRegion;#candidate seed whose middle coordinate is in the "region"
    for my $i (0 .. $#{$read->{QNAME}}){
        my $mid = ($read->{START}[$i]+$read->{END}[$i])/2;#middle coordinate
        if ($mid >= ($bestWin-1)*$opts{w} && $mid <= ($bestWin+1)*$opts{w}){
            push @seedRegion , $i;
        }
    }
    #set region marker
    my %markRegion;#filtered heter-snp-marker in the "region"(sub set of %filter)
    for my $kpos ( keys %filter){
        if ($kpos >= ($bestWin-1)*$opts{w} && $kpos <= ($bestWin+1)*$opts{w}){
            $markRegion{$kpos}++;
        }
    }
    #traverse all reads for alt allel and ref allel
    my $allRead = 0;#all candidate read
    my %snpArray;#all alt-allel snp
    my %iArray;#all alt-allel snp $iArray{i}{pos}
    my %snpFre;#
    my $ind = 0;#index of @seedRegion
    my %index;#value is index of @seedRegion, key is read NO.
    for my $i (@seedRegion){#$i is read NO. ordered in 0..100
        $index{"$i"}="$ind";
        $ind++;
        $allRead++;
        for (my $j=0; $j <= $#{$read->{SNPPOS}[$i]}; $j++){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            if (exists $markRegion{$pos}){
                $snpArray{$pos}{$i} = $alt;
                $snpFre{$pos}{$alt}++;
                $iArray{$i}{$pos} = $alt;
            }
        }
        for my $kpos (keys %markRegion){
            unless(exists $snpArray{$kpos}{$i}){
                $snpArray{$kpos}{$i} = "R";#ref allel but ignore indels
                $snpFre{$kpos}{"R"}++;#ref allel but ignore indels. It is not exactly
                $iArray{$i}{$kpos} = "R";
            }
        }
    }
    #cutoff of SNP
    my %selectPos;#ref:alt allel > cufoff && <(1-cutoff)
    for my $kpos (keys %snpFre){
        my $refFre = 0;
        my $otherFre = 0;
        for my $kalt (keys %{$snpFre{$kpos}}){
            if ($kalt eq "R"){
                $refFre+=$snpFre{$kpos}{$kalt};
            }else{
                $otherFre += $snpFre{$kpos}{$kalt};
            }
        }
        if (($refFre/($refFre+$otherFre)) >= $opts{e} && ($refFre/($refFre+$otherFre)) < (1-$opts{e})){
            $selectPos{$kpos}++;
        }
    }
    my %fre_hash;#the snp pattern of the selected SNPs
    my %fre_hash_i;#the snp pattern of the selected SNPs (value is read NO)
    for my $ki (keys %iArray){
        my $snpPatten;
        for my $kpos (sort {$a<=>$b} keys %{$iArray{$ki}}){
            if (exists $selectPos{$kpos}){
                $snpPatten .= "$iArray{$ki}{$kpos},";
            }
        }
        $fre_hash{$snpPatten}++;
        $fre_hash_i{$snpPatten} .= "$ki,";
    }

    my $maxPatten;
    my $maxFre = 0;
    for my $kpattern (sort keys %fre_hash){
        my $fre = $fre_hash{$kpattern};
        if ($fre > $maxFre){
            $maxPatten = $kpattern;
            $maxFre = $fre;
            my @temp = split /,/,$fre_hash_i{$kpattern};
            for (my $si=0; $si<= $#temp; $si++){
                $$s0[$si] = $temp[$si];
            }
        }
        #print "-- $kpattern : $fre_hash{$kpattern}\n";
        print PT "$kpattern : $fre_hash{$kpattern}\n";
        print LOG "$kpattern : $fre_hash{$kpattern}\n";
    }
    print LOG "Error! Plase set -w larger (now=$opts{w}) and -e more close to 0 (now=$opts{e})" and die if ($maxFre == 0);
    print LOG "Error! Plase set -w larger (now=$opts{w}) and -e more close to 0 (now=$opts{e})" and die if ($maxFre == 1);
    delete $fre_hash{$maxPatten};
    $maxFre = 0;
    for my $kpattern (sort keys %fre_hash){
        my $fre = $fre_hash{$kpattern};
        if ($fre > $maxFre){
            $maxPatten = $kpattern;
            $maxFre = $fre;
            my @temp = split /,/,$fre_hash_i{$kpattern};
            for (my $si=0; $si<= $#temp; $si++){
                $$s1[$si] = $temp[$si];
            }
        }
    }
    print LOG "Error! Plase set -w larger (now=$opts{w}) and -e more close to 0.5 (now=$opts{e})" and die if ($maxFre == 1);
}

sub errorNum{
    my ($i, $x) = @_;
    my %seqError = %$x;
    my $sE = 0;
    for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
        my $pos = $read->{SNPPOS}[$i][$j];
        my $alt = $read->{SNPALT}[$i][$j];
        if(exists $seqError{$pos}{$alt}){
            $sE++;
        }
    }
    return $sE;
}
sub homoNum{
    my ($i, $x) = @_;
    my %homoSnp = %$x;
    my $hS = 0;
    for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
        my $pos = $read->{SNPPOS}[$i][$j];
        my $alt = $read->{SNPALT}[$i][$j];
        if(exists $homoSnp{$pos}{$alt}){
            $hS++;
        }
    }
    return $hS;
}

sub veen {
    my $f1 = "$opts{o}/phase.0.qname";
    my $f2 = "$opts{o}/phase.1.qname";
    my %hash;
    my $p0 = 0;
    my $p1 = 0;
    my $ol = 0;
    open IN , "$f1";
    while(<IN>){
        my ($qname) = (split /\t/)[0];
        $hash{$qname}++;
        $p0++;
    }
    close IN;
    open IN , "$f2";
    while(<IN>){
        my ($qname) = (split /\t/)[0];
        $ol++ if (exists $hash{$qname});
        $p1++
    }
    close IN;
    print "-- Reads of Phase 0: $p0\n";
    print "-- Reads of Phase 1: $p1\n";
    print "-- Same reads between 2 phase:$ol\n";
    print PT "-- Reads of Phase 0: $p0\n";
    print PT "-- Reads of Phase 1: $p1\n";
    print PT "-- Same reads between 2 phase:$ol\n";
}

sub artSeed {
    my ($X, $Y, $Z) = @_;
    my @seed = @$X;
    my $seedSnp = $Y;
    my %filter = %$Z;
    my %tempSnp;
    my $artMark = 0;#markers num in artificial region
    for my $s (@seed){
        for my $j (0 ..$#{$read->{SNPPOS}[$s]}){
            my $pos = $read->{SNPPOS}[$s][$j];
            my $alt = $read->{SNPALT}[$s][$j];
            my $qual = $read->{SNPQUAL}[$s][$j];
            if ($pos > ($bestWin-1)*$opts{w} && $pos <= ($bestWin+1)*$opts{w}){
                $tempSnp{$pos}{$alt} += $qual;
            }
        }
        for my $pos ( (($bestWin-1)*$opts{w}) .. (($bestWin+1)*$opts{w})){
            if (exists $refSnp{$read->{QNAME}[$s]}{$pos}){#this read is ref-allel
                 $tempSnp{$pos}{"R"}+=1-10**(($refSnp{$read->{QNAME}[$s]}{$pos})/(0-10));
            }
        }
    }
    for my $kpos (keys %tempSnp){
        my $maxV = 0;
        my $maxAlt = "NA";
        next if (!exists $filter{$kpos});
        for my $kalt (keys %{$tempSnp{$kpos}}){
            if ($tempSnp{$kpos}{$kalt} > $maxV){
                $maxAlt = $kalt;
                $maxV = $tempSnp{$kpos}{$kalt};
                $artMark++;
            }
        }
        $$seedSnp{$kpos} = $maxAlt;
    }
    return $artMark;
}

########## ########## Help and Information ########## ##########
sub help{
print "*** Phase reads into 2 haplotype, using BAM format files.
Usage:perl thisScript.pl -o $opts{o} -w $opts{w} -u $opts{u} -p $opts{p} -d $opts{d} *.bam
\t-h For more information.
\t-t Delet the TEMP file directory.
\t-c Coincide SNP proportion.(default=$opts{c})
\t-w Window size of seed selection.(default=$opts{w}bp)
\t-o Output directory.(default=$opts{o})
\t-u Phase 0 and 1 > score cutoff. 0.0--1.0(default=$opts{u})
\t-f First seed read NO. 
\t-s Second seed read NO. 
\t-p Upper heter snp cutoff, alt fre/seq depth.(default=$opts{p})
\t-d Lowwer heter snp cutoff, alt fre/seq depth.(default=$opts{d})
\t*.bamFiles: BAM files for small region
";
}

sub info{
print "
*** Describe: This script is designed for 100% covered (well seqenced) region PacBio CCS reads two haplotypes (diploid) phasing. For example, applying HLA/MHC region capture sequencing and then using this script for HLA-A gene region CCS reads phasing. This script performs better than SAMtools phase espcially when 2 haplotypes are in unbalanced sequencing coverage fold.And print out QNAME of each read.Some abnormity long read with little overlap with this region will be discarded.
*** Sugestion(1): Candidate reads pattern:R is ref allel. If the length of pattern is not long enough. You need to set -w larger. If there is no 2 siginificent high value in all value of pattern.
***Sugestion(2): If the \"extend times\" is too less than half of all read number. You need to set seed read NO. by -f -s.
*** Proformance: This script is written in Perl, with lots of C style pointers in functions. With dozen of times optimizing of algorithm, I try to make the phasing result more accurate. One more thing, unfortunately, I used planty of hash data structure to simplify programming, so that the script will not run very fast.
*** Email: zhouze\@genomics.cn
  _________  _________ 
  | __|__ |  |___|___| 
  |___|___|  |___|___| 
  |  ___  |      |     
  | |___| |      |     
  /      \\|      |    
\n";
}
