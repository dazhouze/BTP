#!/usr/bin/perl -w
use strict;
use Getopt::Std;

########## ########## Get paramter ########## ##########
my %opts;
getopts('ho:w:e:c:', \%opts);
$opts{c} = 0.5 unless ($opts{c});#coincide SNP proportion when extending.
$opts{w} = 500 unless ($opts{w});#window size of seed region selection
$opts{e} = 0.3 unless ($opts{e});#cutoff value of seed (SNP) pattern selection

&help and &info and die if ($opts{h}); 
&help and die unless ($opts{o} && @ARGV);

system "mkdir -p $opts{o}";
open LOG , ">$opts{o}/log.txt";#log file
open RES , ">$opts{o}/result.txt";#result file
print RES "NO.\tChr\tStart\tEnd\n";

my $BUFF = 30;
my @buffer;
=from
     $start;
     $end;
     $qname;
     [ @line_snp_pos ];
     [ @line_snp_alt ];
     [ @line_snp_qual ];
=cut
my $if_start = 1;#if start select two seeds 
#best hit read only
open IN , "samtools view $ARGV[0] |" ;

while (<IN>){
    next if (/^@/);
    chomp;
    my $con = $_;
    if ($if_start){
        $if_start = 0;
        my @seedRange = &seedDetermine($con);#every heter snp pos of seed region
        &seedSelect(\@buffer, \@seedRange);  
        &readExtend(\@buffer);  
    }
    die;
    #buffer array read extension
    ########## ########## Read extension ########## ##########
    ########## ########## Read extension ########## ##########
}
my $read;

close IN;
close LOG;
close RES;

########## ########## Functions ########## ##########
sub detectSnp {
    my $con = $_[0];
    chomp($con);
    my %seq_hash;#sequence based (start from 1) hash {ref_pos}{base N.O}{ref}{alt}
    my @line_snp;#all information
    my @line_snp_pos;
    my @line_snp_alt;
    my @line_snp_qual;#qualty of each snp
    
    my ($qname,$chr,$pos,$cigar,$seq,$qual)=(split /\t/,$con)[0,2,3,5,9,10];
    $con =~ /MD:Z:(.*?)\t/;
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
                        push @line_snp_pos , $krp; 
                        push @line_snp_alt , $kb; 
                        push @line_snp_qual , (ord($kqu)-33); 
                    }
                }
            }
        }
    }
    push @line_snp, $start;
    push @line_snp, $end;
    push @line_snp, $qname;
    push @line_snp, [ @line_snp_pos ];
    push @line_snp, [ @line_snp_alt ];
    push @line_snp, [ @line_snp_qual ];
    return @line_snp;
}

sub seedDetermine{
    my ($con) = @_;
    my $lineNum = 0;
    ########## ########## Candidate reads(left most 30 reads) ########## ##########
    for (; $lineNum < $BUFF; $lineNum++){
        my @line_snp = &detectSnp($con);#$start $end $qname @line_snp_pos @line_snp_alt @line_snp_qual
        push @buffer, [ @line_snp ];
        $con = <IN>;
    }
    ########## ########## Seed finding ########## ##########
    my %seedSnpFre;#hash for detect seed snp frequence
    my %seedMark;#heter-snp marker of seed selection
    #alt-allel initialization
    for my $i (0 .. $#buffer){
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            my $alt = $buffer[$i][4][$j]; 
            my $qua = $buffer[$i][5][$j]; 
            $seedSnpFre{$pos}{$alt}++;
        }
    }
    #ref-allel initialization (indels are ignore)
    for my $i (0 .. $#buffer){
        my $start = $buffer[$i][0];
        my $end = $buffer[$i][1];
        #my $qname = $buffer[$i][2];
        my %readSnpPos;
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            $readSnpPos{$pos}++;
        }
        for my $p ($start .. $end){
            if (exists $seedSnpFre{$p} && !exists $readSnpPos{$p}){
                $seedSnpFre{$p}{"R"}++;
            }
        }
    }
    #figure out heter-snp marker
    for my $kpos (sort {$a<=>$b} keys %seedSnpFre){#qname 
        my $allFre = 0;
        for my $kalt (keys %{$seedSnpFre{$kpos}}){#pos 
            $allFre++;
        }
        if( exists $seedSnpFre{$kpos}{"R"} && $seedSnpFre{$kpos}{"R"} >($opts{e}*$allFre) && $seedSnpFre{$kpos}{"R"}<($allFre*(1-$opts{e}))){
            $seedMark{$kpos}++;
            #print "$kpos is heter-snp pos\n";
        }
    }
    #artificial seed
    my @seedRange = reverse sort {$a<=>$b} keys %seedMark;#$seedRange[0] is start and seedRange[1] is end
    return @seeRange;
}


sub seedSelect{
    my ($X, $Y) = @_;#$s0 is \@seed0 $s1 is \@seed1;
    my @buffer = @$X;
    my @seedRange = @$Y;

    #find reads in the "region"
    #set region marker
    my %markRegion;#filtered heter-snp-marker in the "region"(sub set of %filter)
    for my $kpos (@seedRange){
        $markRegion{$kpos}++;
    }
    #traverse all reads for alt allel and ref allel
    my %snpArray;#all alt-allel snp
    my %iArray;#all alt-allel snp $iArray{i}{pos}
    my %snpFre;#
    my $ind = 0;#index of @buffer
    my %index;#value is index of @buffer, key is read NO.
    for my $i (0 .. $#buffer){#$i is read NO. ordered in 0..100
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            my $alt = $buffer[$i][4][$j]; 
            $readSnpPos{$pos}++;
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

sub range_marker{
    my ($read_s , $read_e) = @_;#read start pos and end pos
    my $s = 0;#counter

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

=from

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

    sub markValue{#belong to readExtend part
        my ($i, $X, $Y, $Z, $R, $phaseSnp) = @_;
        my %readSnp0 = %$X;
        my @range = @$R;
        my $numMark = &range_marker(($read->{START}[$i]>$range[0]?$read->{START}[$i]:$range[0]), ($read->{END}[$i]<$range[1]?$read->{END}[$i]:$range[1]));# marker alt should be detect
        my $markYes = 0;#SNP is same to phase_0 marker alt (position and alt-allel)
        my $markNo = 0;#SNP is different to phase_0 marker alt (ref-allel or other base)
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
=cut

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
