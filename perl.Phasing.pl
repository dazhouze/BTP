#!/usr/bin/perl -w
use strict;
use Getopt::Std;

########## ########## Get paramter ########## ##########
my %opts;
getopts('hto:w:m:p:d:s:f:', \%opts);
$opts{m}=0.75 unless ($opts{m});#read overlap with detected region
$opts{u}=0.55 unless ($opts{u});#phase 0 cutoff value of scoring
$opts{w} = 750 unless ($opts{w});#window size = 500 bp
$opts{p} = 0.25 unless ($opts{p});#upper heter snp cutoff, alt fre/seq depth
$opts{d} = 0.75 unless ($opts{d});#lowwer heter snp cutoff, alt fre/seq depth

if ($opts{o}){
    system "mkdir -p $opts{o}";
    system "mkdir -p $opts{o}/TEMP";
}else{
    $opts{o}="./";
    print "Not set the output directory. Auto set:./\n";
}
&help and die if ($opts{h});
die "Phase reads into 2 haplotype, using BAM format files.
Usage:perl thisScript.pl -o ./ -m $opts{m} -w $opts{w} -u $opts{u} -p $opts{p} -d $opts{d} *.bam
\t-h For more information.
\t-t Delet the TEMP file directory.
\t-w Window size of seed selection.(default=$opts{w}bp)
\t-o Output directory.(default=$opts{o})
\t-m min-overlap between read and detected region: 0.0--1.0.(default=$opts{m})
\t-u phase_0 > score cutoff: 0.0--1.0 .(default=$opts{u})
\t-f First seed read NO. 
\t-s Second seed read NO. 
\t-p Upper heter snp cutoff, alt fre/seq depth.(default=$opts{p})
\t-d Lowwer heter snp cutoff, alt fre/seq depth.(default=$opts{d})
\t*.bamFiles: BAM files for small region\n" unless ($opts{o} && $opts{m} && @ARGV);

open OUT , ">$opts{o}/TEMP/score.temp";
open HIT , ">$opts{o}/TEMP/hit.temp";#read extend marker hit
open MK , ">$opts{o}/TEMP/markerSNP.temp";#marker SNP
open PT , ">$opts{o}/TEMP/seedSelect.temp";
#open PH , ">$opts{o}/TEMP/phaseSNP.temp";#phase 0 SNP

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

########## ########## Detect all SNP in all read ########## ##########
my @bam = @ARGV;
&detectSnp("mismatch");
print "   Finish (1/14). BAM file SNP detection.\n";

########## ########## Identify heterozygous SNP marker, seq error and homo SNP ########## ##########

my %filter;#filtered marker hash
my %seqError;#read snp that seem to be sequencing error
my %homoSnp;#read snp that seem to be homozygous snp

my $heterSnpNum = &traverseSnp(\%filter, \%seqError, \%homoSnp);
print "   Finish (2/14). Heterozygous SNP markers:$heterSnpNum determination.\n";

########## ########## Detect ref-allel SNP in all read ########## ##########
my %refSnp;
&detectSnp("match");
print "   Finish (3/14). Ref-allel SNP detection.\n";

########## ########## Set seed ########## ##########
my $bestWin = &traverseForSeedRegion(\%filter);
print "   Finish (4/14). Seed windows region setting. Window:$bestWin.\n";
print "-- Seed selection hash.\n-- Key : Value\n";

my ($seed0,$seed1) = &seedSelect($bestWin, \%filter, \%seqError, \%homoSnp);#return 2 shortest read in diff hap as seed.
print "   Finish (5/14). Auto seed setting. Phase 0 seed NO: $seed0 length: $read->{LEN}[$seed0]. Phase 1 seed NO: $seed1 length: $read->{LEN}[$seed1].\n";

$seed0 = $opts{f} if ($opts{f});
$seed1 = $opts{s} if ($opts{s});

&Phase($seed0,"phase.0",5,0);
&Phase($seed1,"phase.1",5,1);

sub Phase{
    my ($seed, $fileName,$step,$pha) = @_;
    ########## ########## Initialize phase_0 SNP ########## ##########
    my %seedSnp;# snp of seed read
    my %phase0Snp;
    my $allRead = $#{$read->{LEN}} + 1;

    &phase0initial(\%phase0Snp, \%filter, \%refSnp);#initial heter snp marker pos to phase0
    &seedInitial(\%phase0Snp, \%filter, \%refSnp, $seed);#initial seed to phase 0
    $step--;
    print "   Finish ($step/14). Phase $pha SNPs initialization.\n";

    ########## ########## Grow SNP tree (phase 0 SNP markers) ########## ##########
    my @range;#detect range of genome [0]:start pos, [1]:end pos
    $range[0] = $read->{START}[$seed];
    $range[1] = $read->{END}[$seed];
    my $extendTime = &readExtend($seed, \@range, \%phase0Snp, \%filter, \%refSnp, \%seqError, \%homoSnp);
    $step++;
    print "-- Phase $pha seed: read NO.$seed extend times: $extendTime (all read: $allRead)\n";
    print "   Finish ($step/14). Phase $pha heter-SNP-marker tree growth.\n";

    ########## ########## Filter phase_0 heter SNP markers ########## ##########
    my %phase0SnpFilter;
    &markerFilter(\%phase0Snp, \%phase0SnpFilter);
    $step++;
    print "   Finish ($step/14). Phase $pha heter SNP markers filtering.\n";

    ########## ########## Scoring all reads ########## ##########
    my %qnameMark;#hash of qname and markYes(hit marker) value
    my @markerVal;#arrary of marker Yes(hit) value for midium caculation
    &scoring(\%phase0SnpFilter, \@markerVal, \%qnameMark, \%filter, \%seqError, \%homoSnp, \%refSnp);
    $step++;
    print "   Finish ($step/14). All reads scoring.\n";

    ########## ########## Determine 2 haplotype ########## ##########
    &printResult(\@markerVal, \%qnameMark, "$fileName.qname");
    $step++;
    print "   Finish ($step/14). Phase $pha result printing.\n";
}

close OUT;
#close PH;
close MK;
close HIT;
system"rm -r $opts{o}/TEMP/" if ($opts{t});
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
            die "No MD field was found in SAM/BAM file." unless($md);
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
                    die "Unknow Cigar character.\n";
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
                                        print"Script error:$ksp $krp $kt == M\n" if($kt eq "I");
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
                                print"Script error:$krp\t$ksp\t$kb\t$kt\t$kr\n" if ($kt eq "Mis" && $kb eq $kr);
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
    my $phase0Snp = shift;
    my $kpos = shift;

    my $maxAlt = "NA";
    my $max = 1;
    for my $kalt (keys %{$$phase0Snp{$kpos}}){
        if ($$phase0Snp{$kpos}{$kalt} > $max){
        $max = $$phase0Snp{$kpos}{$kalt};
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
    my ($x, $y, $z, $seed) = @_;
    my $phase0Snp = $x;
    my %filter = %$y; 
    my %refSnp = %$z;
    for my $j (0 ..$#{$read->{SNPPOS}[$seed]}){
        my $pos = $read->{SNPPOS}[$seed][$j];
        my $alt = $read->{SNPALT}[$seed][$j];
        my $qual = $read->{SNPQUAL}[$seed][$j];
        if (exists $filter{$pos}){#same to filtered marker
            $$phase0Snp{$pos}{$alt} += $qual;
        }
    }
    for my $z ($read->{START}[$seed] .. $read->{END}[$seed]){
        if (exists $refSnp{$z}){
            $$phase0Snp{$z}{"ref"}+=1;
        }
    }
}

sub phase0initial{
    my ($x,$y,$z)= @_;
    my $phase0Snp = $x;
    my %filter = %$y;
    my %refSnp = %$z;
    for my $i (0 .. $#{$read->{QNAME}}){#initial phase 0 SNP of seed read
        for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            my $qual = $read->{SNPQUAL}[$i][$j];
            if (exists $filter{$pos}){#same to filtered marker
                $$phase0Snp{$pos}{$alt}+=(1/300);
            }
        }
    }
    for my $kpos (sort {$a<=>$b} keys %filter){
        if (exists $refSnp{$kpos}){
            $$phase0Snp{$kpos}{"ref"}+=(1/300);
        }
    }
}

sub readExtend{
    my ($a, $b, $c, $d, $e, $f, $g) = @_;
    my $seed = $a;
    my @range = @$b;
    my $phase0Snp = $c;
    my %filter = %$d;
    my %refSnp = %$e;
    my %seqError = %$f;
    my %homoSnp = %$g;

    my $maxredo = $#{$read->{QNAME}};#read number
    my $p0r = 0;#read number determined as phase_0
    my %readConsider;# read been considered
    $readConsider{$seed}++;# read been considered
    while ($maxredo){
        $maxredo--;
        for my $i (0 .. $#{$read->{QNAME}}){
            unless(exists $readConsider{$i}){ # skip phased read
                if (&range_overlap($range[0],$range[1],$read->{START}[$i],$read->{END}[$i])<=($opts{m}*$read->{LEN}[$i])){# skip no overlap read(overlap < 0.4 length of read)and check it later (redo)
                    next;
                }else{
                    $readConsider{$i}++;
                    my $ft_re_num = 0;#filtered SNP marker number overlap with read snp
                    my $numMark = &range_marker(($read->{START}[$i]>$range[0]?$read->{START}[$i]:$range[0]), ($read->{END}[$i]<$range[1]?$read->{END}[$i]:$range[1]));# marker alt should be detect
                    my $markYes = 0;#SNP is same to phase_0 marker alt (position and alt-allel)
                    my $markNo = 0;#SNP is different to phase_0 marker alt (ref-allel or other base)
                    my $sE = 0; #sequencing error
                    my $hS = 0; #homo snp 
                    my $otherSnp = 0;#SNP is not exist in marker alt
                    my $rS = 0;#SNP is heterzogous one is snp one is ref
                    my %readSnp;
                    my %readSnp0;
                    for my $j (0 ..$#{$read->{SNPPOS}[$i]}){
                        my $pos = $read->{SNPPOS}[$i][$j];
                        my $alt = $read->{SNPALT}[$i][$j];
                        my $qual = $read->{SNPQUAL}[$i][$j];
                        $readSnp{$pos}{$alt}=$qual;
                        $readSnp0{$pos}=$alt;
                        if(exists $filter{$pos}){
                            if ($alt eq &max_alt($phase0Snp, $pos)){#pos and alt is same(alt-allel)
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
                    for my $pos (($read->{START}[$i]>$range[0]?$read->{START}[$i]:$range[0]) .. ($read->{END}[$i]<$range[1]?$read->{END}[$i]:$range[1])){
                    #for my $pos ($read->{START}[$i] .. $read->{END}[$i]){
                        if (exists $refSnp{$pos}){
                            my $alt = "ref";
                            if ($alt eq &max_alt($phase0Snp, $pos)){#pos and alt is same(alt-allel)
                                $markYes++;
                            }
                            else{#ref-allel
                                $markNo++;
                            }
                        }
                    }
                    #judgement of phase_0/phase_1
                    my $xAxis=-1;
                    my $yAxis=-1;
                    #print "other phase seed index $i $read->{LEN}[$i]\n";
                    if ($markYes+$markNo+$sE){
                        $xAxis=$markYes/($markYes+$markNo+$sE+$hS);
                        $yAxis=$markNo/($markYes+$markNo+$sE+$hS);
                    }
                    #if($markYes >= $markNo){#phase0
                    if(0.4*$xAxis >= 0.6*$yAxis){#phase0
                        $p0r++;
                        $range[0]=($range[0]<$read->{START}[$i]?$range[0]:$read->{START}[$i]);#re-new start of range
                        $range[1]=($range[1]>$read->{END}[$i]?$range[1]:$read->{END}[$i]);#re-new end of range
                        for my $z ($read->{START}[$i] .. $read->{END}[$i]){#z is position from read start to end
                            if (exists $filter{$z}){#there should be SNP marker
                                if (exists $readSnp{$z}){#there is read SNP
                                    for my $kalt (keys %{$readSnp{$z}}){
                                        #SNP base qual -10log10{Pr}
                                        my $kqual = 1-10**($readSnp{$z}{$kalt}/(0-10));#base qual of SNP (<1)
                                        $$phase0Snp{$z}{$kalt}+=$kqual;
                                    }
                                }
                                elsif (exists $refSnp{$read->{QNAME}[$i]}{$z}){#this read is ref-allel
                                    $$phase0Snp{$z}{"ref"}+=1-10**(($refSnp{$read->{QNAME}[$i]}{$z})/(0-10));
                                }
                                else{#complex situation
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    undef %readConsider;
    return $p0r;
}

sub markerFilter{
    my ($x, $y) = @_;
    my %phase0Snp = %$x;
    my $phase0SnpFilter = $y;
    for my $kpos (keys %phase0Snp){
        my $max = 0.1;
        my $maxAlt = "NA";
        for my $kalt (sort keys %{$phase0Snp{$kpos}}){
            if ($phase0Snp{$kpos}{$kalt} > $max){
                $max = $phase0Snp{$kpos}{$kalt};
                $maxAlt = $kalt;
            }
        }
        $$phase0SnpFilter{$kpos}="$maxAlt";
        #print PH "$kpos\t$maxAlt\n";
    }
}

sub scoring{
    my ($a, $b, $c, $d, $e, $f, $g) = @_;
    my %phase0SnpFilter = %$a;
    my $markerVal = $b;
    my $qnameMark = $c;
    my %filter = %$d;
    my %seqError = %$e;
    my %homoSnp = %$f;
    my %refSnp = %$g;

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
                    if ($alt eq $phase0SnpFilter{$z}){#read snp alt is same to phase0 library
                        $markYes+=$kqual;
                    }
                    elsif(exists $seqError{$z}{$alt}){#read snp is sequencing error
                        $sE++;
                    }
                    elsif(exists $homoSnp{$z}{$alt}){#read snp is homo-snp which is not marker
                        $hS++;
                    }
                }
                elsif($phase0SnpFilter{$z} eq "ref"){#there is no read SNP because SNP marker is ref-alt
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
        print OUT "$read->{LEN}[$i]\t$numMark\t$markYes\t$markNo\t$sE\t$hS\n";
        if ($numMark){
            $markYes=$markYes/$numMark; 
            $markNo=$markNo/$numMark;
        }else{
            $markYes=1;#if there is no heter SNP, which means can not judge, consider this read belongs to all 2 hap 
        }
        push @$markerVal , $markYes;
        $$qnameMark{$read->{QNAME}[$i]}=$markYes;
    }
}

sub printResult{
    my ($a, $b, $fileName) = @_;
    my @markerVal = @$a;
    my %qnameMark = %$b;

    my $count = @markerVal;
    my $p0=0;#reads number of phase_0
    my $p1=0;#reads number of phase_1
    #my $markerCutOff0 = $sortMarkerVal[0]*0.45;#phase0 > 
    #my $markerCutOff1 = $sortMarkerVal[0]*0.20;#phase0 > 
    #my $markerCutOff1 = $sortMarkerVal[0]*0.20;#phase1 <
    #my @reSortMarkerVal = reverse sort @markerVal;#from small to large
    my @sortMarkerVal =  sort @markerVal;#from small to large
    my $markerCutOff0 = $sortMarkerVal[int($opts{u}*($count))];#small
    print "-- Marker hit cutoff value: $markerCutOff0\n";
    #####;
    open OUT0 , ">$opts{o}/$fileName";
    for my $kq (keys %qnameMark){
        my $markYes = $qnameMark{$kq};
        if ($markYes >= $markerCutOff0){#phase0
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
    my ($bestWin, $x, $y, $z) = @_;
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
    my %markRegion;#filtered heter-snp-marker in the "region"
    for my $kpos ( keys %filter){
        if ($kpos >= ($bestWin-1)*$opts{w} && $kpos <= ($bestWin+1)*$opts{w}){
            $markRegion{$kpos}++;
        }
    }
    #traverse all reads for alt allel and ref allel
    my $allRead = 0;#all candidate read
    my %snpArray;
    my %halfSnp;
    my $ind = 0;#index of @seedRegion
    my %index;#value is index of @seedRegion, key is read NO.
    my %fre_hash;#the snp patten of the selected SNPs
    my %fre_hash_i;#the snp patten of the selected SNPs
    for my $i (@seedRegion){#$i is read NO. ordered in 0..100
        $index{"$i"}="$ind";
        $ind++;
        $allRead++;
        for (my $j=0; $j <= $#{$read->{SNPPOS}[$i]}; $j++){
            my $pos = $read->{SNPPOS}[$i][$j];
            my $alt = $read->{SNPALT}[$i][$j];
            if (exists $markRegion{$pos}){
                $snpArray{$pos}{$i} = $alt;
                $halfSnp{$pos}{$alt}++;
            }
        }
        for my $kpos (keys %markRegion){
            unless(exists $snpArray{$kpos}{$i}){
                $snpArray{$kpos}{$i} = "R";#ref allel
                $halfSnp{$kpos}{"R"}++;
            }
        }
        my $snpPatten;
        for my $kpos (sort {$a<=>$b} keys %snpArray){
            $snpPatten .= "$snpArray{$kpos}{$i},";
        }
        $fre_hash{$snpPatten}++;
        $fre_hash_i{$snpPatten} .= "$i,";
    }
    my $maxPatten;
    my $maxFre = 0;
    my @clu0;
    my @clu1;
    for my $kpatten (sort keys %fre_hash){
        my $fre = $fre_hash{$kpatten};
        if ($fre > $maxFre){
            $maxPatten = $kpatten;
            $maxFre = $fre;
            @clu0 = split /,/,$fre_hash_i{$kpatten};
        }
        print "-- $kpatten : $fre_hash{$kpatten}\n";
        print PT "$kpatten : $fre_hash{$kpatten}\n";
    }
    delete $fre_hash{$maxPatten};
    print "-- Phase 0 candidate reads patten: $maxPatten NO:$fre_hash_i{$maxPatten}\n";
    $maxFre = 0;
    for my $kpatten (sort keys %fre_hash){
        my $fre = $fre_hash{$kpatten};
        if ($fre > $maxFre){
            $maxPatten = $kpatten;
            $maxFre = $fre;
            @clu1 = split /,/,$fre_hash_i{$kpatten};
        }
    }
    print "-- Phase 1 candidate reads patten: $maxPatten NO:$fre_hash_i{$maxPatten}\n";

    my $s = 0;
    my $minLen;
    my $seed0;
    my $seed1;
    for my $kci (@clu0){
        if($s == 0){
            $seed0 = $kci;
            $minLen = $read->{LEN}[$kci];
        }elsif($read->{LEN}[$kci] < $minLen){
            $seed0 = $kci;
            $minLen = $read->{LEN}[$kci];
        }
        $s++;
        #my $sE=&errorNum($kci, \%seqError);
        #my $hS=&homoNum($kci, \%homoSnp);
        #print "$s\t$kci\t$read->{LEN}[$kci]\t$sE\t$hS\n";
    }
    $s = 0;
    for my $kci (@clu1){
        if($s == 0){
            $seed1 = $kci;
            $minLen = $read->{LEN}[$kci];
        }elsif($read->{LEN}[$kci] < $minLen){
            $seed1 = $kci;
            $minLen = $read->{LEN}[$kci];
        }
        $s++;
        #print "$s\t$kci\tclu1\n";
    }
    return $seed0, $seed1;
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

########## ########## End ########## ##########

sub help{
print "
** Describe: This script is designed for 100% covered (well seqenced) region PacBio CCS reads two haplotypes (diploid) phasing. For example, applying HLA/MHC region capture sequencing and then using this script for HLA-A gene region CCS reads phasing. This script performs better than SAMtools phase espcially when 2 haplotypes are in unbalanced sequencing coverage fold.And print out QNAME of each read.Some abnormity long read with little overlap with this region will be discarded.
** Sugestion(1): Candidate reads patten:R is ref allel. If the length of patten is not long enough. You need to set -w larger.
** Sugestion(2): If the \"extend times\" is too less than half of all read number. You need to set seed read NO. by -f -s.
** Proformance: This script is written in Perl, with lots of C style pointers in functions. With dozen of times optimizing of algorithm, I try to make the phasing result more accurate. One more thing, unfortunately, I used planty of hash data structure to simplify programming, so that the script will not run very fast.
** Email: zhouze\@genomics.cn
  _________  _________ 
  | __|__ |  |___|___| 
  |___|___|  |___|___| 
  |  ___  |      |     
  | |___| |      |     
  /      \\|      |    
\n";
}
