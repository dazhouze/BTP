#!/usr/bin/perl -w
use strict;
use Getopt::Std;

########## ########## Get paramter ########## ##########
my %opts;
getopts('ho:w:e:c:d:', \%opts);
$opts{c} = 0.95 unless ($opts{c});#coincide SNP proportion when extending.
#$opts{w} = 500 unless ($opts{w});#window size of seed region selection
$opts{e} = 0.45 unless ($opts{e});#cutoff value of seed (SNP) pattern selection
$opts{d} = 0.6 unless ($opts{d});#depth cutoff value of seed position selection

&help and &info and die if ($opts{h}); 
&help and die unless ($opts{o} && @ARGV);

system "mkdir -p $opts{o}";
open LOG , ">$opts{o}/log.txt";#log file
open RES , ">$opts{o}/result.txt";#result file
print RES "NO.\tChr\tStart\tEnd\n";

my $BUFF = 100;
my @buffer;#stack

#data struct
#$start;
#$end;
#$qname;
#[ @line_snp_pos ];
#[ @line_snp_alt ];
#[ @line_snp_qual ];

my @phase0;#reads NO.
my @phase1;#reads NO.
my @range0;#detect range of phase0
my @range1;#detect range of phase1

my $read_status = 1;#1:need to select two seeds. 0:no need to select seeds. 2:undeiside 
#best hit read only
open FH , "samtools view $ARGV[0] |" ;
my $l = 0;
while (1) {
    ##my @lines ;#content of each line

    if ($read_status){# begin of a porces local *FH file handle willnot change it

        $read_status = 0;
        ########## ########## put candidate reads in buffer ########## ##########
        for (my $lineNum = 0; $lineNum < $BUFF; $lineNum++){
            my $con = &read_text(*FH);
            if ($con eq "Not_MHC"){
                $lineNum--; 
                next;
            }
            my @line_snp = &detectSnp($con);#$start $end $qname @line_snp_pos @line_snp_alt @line_snp_qual
            push @buffer, [ @line_snp ];
        }
        #print "buffer: $buffer[20][2]\n";

        ########## ########## figure out heter-snp pos of artifical seed ########## ##########
        my @seedRange;
        &seedPosDetermine(\@buffer, \@seedRange);#every heter snp pos of seed region
        #print "select pos:@seedRange\n";
        print "seed range: ", $seedRange[$#seedRange]-$seedRange[0], "bp\n";

        ########## ########## candidate seed reads selection ########## ##########
        my @seed0;#read NO. array
        my @seed1;#read NO. array
        &readDetermine(\@buffer, \@seedRange, \@seed0, \@seed1); 
        print "seed 0 :@seed0\nseed 1:@seed1\n";

        ########## ########## construct artifical seed ########## ##########
        my @seed0Snp;#artifical seed read and snp [pos][alt]
        my @seed1Snp;#artifical seed read and snp [pos][alt]
        &artificialSeed(\@seed0, \@seed0Snp, \@seedRange, "s0"); 
        &artificialSeed(\@seed1, \@seed1Snp, \@seedRange, "s1"); 
        #print "seed 0 :@seed0Snp\nseed 1:@seed1Snp\n";

        ########## ########## Read extension (extend the first buffer arrya)########## ##########
        &readExtend(\@seed0Snp, \@phase0, \@buffer, \@range0);  
        &readExtend(\@seed0Snp, \@phase1, \@buffer, \@range1);  
        die;
    }
    #else {#
    #    ########## ########## Read extension ########## ##########
    #    #buffer array read extension
    #    &dataAccess((*FH), \@buffer);
    #    my $con = &read_text(*FH);
    #    $l++;
    #    #print "$qname\n";
    #    die if ($l>10);
    #    #if(){
    #    #   $read_status = 1;#meet break point of a hap
    #    #   next;
    #    #}
    #}

    #loop end check
    last if  (&ifEOF(*FH));
}
my $read;

close FH;
close LOG;
close RES;

########## ########## Functions ########## ##########

sub max_alt{
    my ($X, $Y) = @_;
    my $hash = $X;
    my $key = $Y; 

    my $maxAlt = "0";
    my $max = 0;
    for my $kalt (keys %{$$hash{$key}}){
        if ($$hash{$key}{$kalt} > $max){
        $max = $$hash{$key}{$kalt};
        $maxAlt = $kalt;
        }
    }
    return $maxAlt;
}

sub artificialSeed{
    my ($X, $Y, $Z, $W) = @_;#\@seed0, \@seed1Snp \@seedRange $seed_name
    my @seed = @$X;
    #my $seed_name = $W;
    my %pos_hash;
    my @seedRange = @$Z;
    for my $pos (@$Z){#@seedRange = @$Z;#positions of SNPs
        $pos_hash{$pos}++;
    }
    my %allReadSnp;
    for my $i (0 .. $#seed){#$i is read NO. ordered in 0..100
        my %readSnp;#singel read SNP
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            my $alt = $buffer[$i][4][$j]; 
            $readSnp{$pos}{$alt}++;
        }
        for my $kpos (keys %pos_hash){
            if(!exists $readSnp{$kpos}){
                $readSnp{$kpos}{"R"}++;
            }
        } 
        for my $kpos (keys %readSnp){
            for my $kalt (keys %{$readSnp{$kpos}}){
                $allReadSnp{$kpos}{$kalt}++;
            }
        }
    }

    $$Y[0] = $seedRange[0]; 
    $$Y[1] = $seedRange[($#seedRange)]; 
    $$Y[2] = $W; 
    my @line_snp_pos;
    my @line_snp_alt;
    my @line_snp_qual;
    for my $pos (@$Z){#@seedRange
        push @line_snp_pos, $pos;
        push @line_snp_alt, &max_alt(\%allReadSnp, $pos);
        push @line_snp_qual, 1;
    }
    $$Y[3] = [@line_snp_pos]; 
    $$Y[4] = [@line_snp_alt]; 
    $$Y[5] = [@line_snp_qual]; 
}

sub readDetermine{
    my ($X, $Y, $s0, $s1) = @_;#$s0 is \@seed0 $s1 is \@seed1;
    my @buffer = @$X;
    my @seedRange = @$Y;

    #find reads in the "region"
    #set region marker
    my %selectPos;#filtered heter-snp-marker in the "region"(sub set of %filter)
    for my $kpos (@seedRange){
        $selectPos{$kpos}++;
    }
    my %fre_hash;
    my %fre_hash_i;
    for my $i (0 .. $#buffer){#$i is read NO. ordered in 0..100
        #$start = $buffer[$i][0][$j]; $end = $buffer[$i][1][$j];
        if ($seedRange[0] >= $buffer[$i][0] && $seedRange[1] <= $buffer[$i][1]){#only consider seed range covered reads 
            my $snpPatten;
            my %readSnp;
            for my $j (0 .. $#{$buffer[$i][3]}){
                my $pos = $buffer[$i][3][$j]; 
                my $alt = $buffer[$i][4][$j]; 
                my $qua = $buffer[$i][5][$j]; 
                $readSnp{$pos}=$alt;
            }
            for my $kpos (@seedRange){#seedRange 
                #print "pos:$pos alt:$alt\n";
                if (exists $readSnp{$kpos}){
                    $snpPatten .= "$readSnp{$kpos},";
                }
                else{#ref-allel
                    $snpPatten .= "R,";
                }
            }
            $fre_hash{$snpPatten}++;
            $fre_hash_i{$snpPatten} .= "$i,";
        }
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
    #print LOG "Error! Plase set -w larger (now=$opts{w}) and -d more close to d (now=$opts{d})" and die if ($maxFre == 0);
    #print LOG "Error! Plase set -w larger (now=$opts{w}) and -d more close to d (now=$opts{d})" and die if ($maxFre == 1);
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
    print "Error! Plase set -w larger (now=$opts{w}) and -d more close to 0.5 (now=$opts{d})" and die if ($maxFre == 1);
    print LOG "Error! Plase set -w larger (now=$opts{w}) and -d more close to 0.5 (now=$opts{d})" and die if ($maxFre == 1);
}

sub seedPosDetermine{#detect apropriate seed snps pos
    my ($X, $Y) = @_;
    my @buffer = @$X;
    #my @seedRange = @$Y;

    my %seedSnpFre;#hash for detect seed snp frequence
    my %seedMark;#heter-snp marker of seed selection
    my %seqDepth;#sequencing depth of each window
    #alt-allel initialization
    for my $i (0 .. $#buffer){
        ########## ########## Snp finding ########## ##########
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            my $alt = $buffer[$i][4][$j]; 
            my $qua = $buffer[$i][5][$j]; 
            $seedSnpFre{$pos}{$alt}++;
        }
    }
    #for my $kw (sort {$a<=>$b} keys %seqDepth){
    #    #print "$kw,$seqDepth{$kw}\n";
    #}
    #ref-allel initialization (indels are ignore)
    for my $i (0 .. $#buffer){
        my $start = $buffer[$i][0];
        my $end = $buffer[$i][1];
        #my $qname = $buffer[$i][2];
        #print "$i,$buffer[$i][0],$buffer[$i][1]\n";
        my %readSnpPos;
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$i][3][$j]; 
            $readSnpPos{$pos}++;
        }
        for my $p ($start .. $end){
            if (exists $seedSnpFre{$p} && !exists $readSnpPos{$p}){# !exists $readSnpPos{$p} it is not a alt-allel but need to add to all hash to determin frequency
                $seedSnpFre{$p}{"R"}++;
            }
        }
    }
    #judice #figure out heter-snp marker
    for my $kpos (sort {$a<=>$b} keys %seedSnpFre){#qname 
        my $allFre = 0;
        for my $kalt (keys %{$seedSnpFre{$kpos}}){#pos 
            $allFre += $seedSnpFre{$kpos}{$kalt};
        }

        next if( $allFre <= ($opts{d}*($BUFF)));
        #print "$kpos $allFre\n";

        if( exists $seedSnpFre{$kpos}{"R"} && $seedSnpFre{$kpos}{"R"} >= ($opts{e}*$allFre) && $seedSnpFre{$kpos}{"R"} <= ($allFre*(1-$opts{e}))){
            $seedMark{$kpos}++;
            #print "$kpos is heter-snp pos\n";
        }
    }
    #add pos in seedRange array
    my $y = 0;#index of @seedRengion
    for my $kpos (sort {$a<=>$b} keys %seedMark){
        $$Y[$y]=$kpos; 
        $y++;
    }
    return 0;
}

sub ifMHC {
    my $con = shift;
    my ($qname,$chr,$pos,$cigar,$seq,$qual)=(split /\t/, $con)[0,2,3,5,9,10];
    my $len = length $seq;
    if ($chr eq "chr6"){
        unless ($pos>33448354 ){
            unless (($pos+$len)<28477797){
                return 1;
            }
        }
    }
    return 0;
}

sub ifEOF{
    local (*FH) = shift;
    return 1 if (eof(*FH));
    return 0;
}

sub read_text{
    local (*FH) = shift;
    my $con = <FH>;
    chomp($con);
    if (&ifMHC($con)) {#only covered MHC region
        #@$lines = split /\t/, $con;
        return $con;
    }
    else{
        return "Not_MHC";
    }   
}

sub detectSnp {
    my $con = $_[0];
    #chomp($con);
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
                $seq_hash{$now_pos}{$now_base}{"D"}{$bQua}{"D"}="delet";
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

sub range_marker{
    #my $numMark = &range_marker(($buffer[$i][0]>$range[0]?$buffer[$i][0]:$range[0]), ($buffer[$i][1]<$range[1]?$buffer[$i][1]:$range[1]));# marker alt should be detect
    my ($read_s , $read_e) = @_;#read start pos and end pos
    my $s = 0;#counter

    return $s;
}

sub readExtend{
    my ($X, $Y, $Z, $W) = @_;#\@seed0Snp, @seed1Snp, \@buffer
    my @range = @$W;

    my $p0r = 0;#read number determined as phase_0
    my $pre_p0r = -1;#previous p0r
    my %readConsider;# read been considered
    my $ovl = 1; 
    my $rC = 1;#read been considered number
    for (my $maxRedo=1+$#buffer; $maxRedo; $maxRedo--){
        my %accordantRead;#reads can be add to "detect region" {consistency}{i}
        last if ($ovl < 0.1);#the overlap is too low
        $rC = keys %readConsider;
        last if ($rC == (1+$#buffer));
        print "range(bp)",$range[1]-$range[0],"\tread consider:$rC\tread extend:$p0r","\toverlap",$ovl,"\n";
        for my $i (0 .. $#buffer){
            if (exists $readConsider{$i}) { # skip phased read
                next;
            }
            else{
                next if (&range_overlap($range[0],$range[1],$buffer[$i][0],$buffer[$i][1])<=($ovl*($buffer[$i][1]-$buffer[$i][0])));
                $readConsider{$i}++;#
                my %conReadSnp0;
                for my $j (0 .. $#{$buffer[$i][3]}){
                    my $pos = $buffer[$i][3][$j]; 
                    my $alt = $buffer[$i][4][$j]; 
                    my $qua = $buffer[$i][5][$j]; 
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
        my $numMark = &range_marker(($buffer[$i][0]>$range[0]?$buffer[$i][0]:$range[0]), ($buffer[$i][1]<$range[1]?$buffer[$i][1]:$range[1]));# marker alt should be detect
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
        #ref-allel ???????????????????????????????
        #for my $pos (($buffer[$i][0]>$range[0]?$buffer[$i][0]:$range[0]) .. ($buffer[$i][1]<$range[1]?$buffer[$i][1]:$range[1])){
        #    if (exists $refSnp{$read->{QNAME}[$i]}{$pos}){
        #        my $alt = "R";
        #        if ($alt eq &max_alt($phaseSnp, $pos)){#pos and alt is same(alt-allel)
        #            $markYes++;
        #        }
        #        else{#ref-allel
        #            $markNo++;
        #        }
        #    }
        #}
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
        
        $$range[0]=($$range[0]<$buffer[$ki][0]?$$range[0]:$buffer[$ki][0]);#re-new start of range
        $$range[1]=($$range[1]>$buffer[$ki][1]?$$range[1]:$buffer[$ki][1]);#re-new end of range
        #print HIT "$pha\t$ki\t$read->{LEN}[$ki]\t$numMark\t$markYes\t$markNo\t$sE\t$hS\n";
        my %readSnp;
        my %readSnp0;
        for my $j (0 .. $#{$buffer[$i][3]}){
            my $pos = $buffer[$ki][3][$j]; 
            my $alt = $buffer[$ki][4][$j]; 
            my $qua = $buffer[$ki][5][$j]; 
            $readSnp{$pos}{$alt}=$qual;
            $readSnp0{$pos}=$alt;
        }
        for my $z ($buffer[$ki][0] .. $buffer[$ki][1]){#z is position from read start to end
            if (exists $filter{$z}){#there should be SNP marker
                if (exists $readSnp{$z}){#there is read SNP
                    for my $kalt (keys %{$readSnp{$z}}){
                        #SNP base qual -10log10{Pr}
                        my $kqual = 1-10**($readSnp{$z}{$kalt}/(0-10));#base qual of SNP (<1)
                        $$phaseSnp{$z}{$kalt}+=$kqual;
                    }
                }
                # ref-allel ?????????????????????????
                #elsif (exists $refSnp{$read->{QNAME}[$ki]}{$z}){#this read is ref-allel
                #    $$phaseSnp{$z}{"R"}+=1-10**(($refSnp{$read->{QNAME}[$ki]}{$z})/(0-10));
                #}
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
    for my $j (0 .. $#{$buffer[$i][3]}){
        my $pos = $buffer[$i][3][$j]; 
        my $alt = $buffer[$i][4][$j]; 
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
    for my $j (0 .. $#{$buffer[$i][3]}){
        my $pos = $buffer[$i][3][$j]; 
        my $alt = $buffer[$i][4][$j]; 
        if(exists $homoSnp{$pos}{$alt}){
            $hS++;
        }
    }
    return $hS;
}

########## ########## Help and Information ########## ##########
sub help{
print "*** Phase reads into 2 haplotype, using BAM format files.
Usage:perl thisScript.pl -o $opts{o} -w $opts{w} *.bam
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
