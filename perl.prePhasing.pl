#!/usr/bin/perl -w
use strict;
use Getopt::Std;

########## ########## Get paramter ########## ##########
my %opts;
getopts('hto:w:p:d:s:f:c:e:u:', \%opts);
$opts{c} = 0.5 unless ($opts{c});#coincide SNP proportion when extending.
$opts{w} = 500 unless ($opts{w});#window size of seed region selection
$opts{e} = 0.3 unless ($opts{e});#cutoff value of seed (SNP) pattern selection

&help and &info and die if ($opts{h}); 
&help and die unless ($opts{o} && @ARGV);

open LOG , ">log.txt";#log file
open RES , ">result.$ARGV[0].txt";#result file
print RES "NO.\tChr\tStart\tEnd\n";

close LOG;
close RES;


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
        print "range(bp)",$range[1]-$range[0],"\tread consider:$rC\tread extend:$p0r","\toverlap",$ovl,"\n";
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
\t-w Window size of seed selection.(default=$opts{w}bp)
\t-o Output directory.(default=$opts{o})
\tbamFile: sorted best hit only BAM file.
";
}

sub info{
print "
*** Describe: This script is designed to figure out the coordinates of regions where two haplotype sequencing depth is equal and good enough (>3X). 
*** Email: zhouze\@genomics.cn
  _________  _________ 
  | __|__ |  |___|___| 
  |___|___|  |___|___| 
  |  ___  |      |     
  | |___| |      |     
  /      \\|      |    
\n";
}
