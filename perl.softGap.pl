#!/usr/bin/perl -w
use strict;

#base (ATGC) in MD field is reference

my $o = "mhc_flank.depth.txt";

open OUT , ">$o" and print "Output file: $o\n";
print OUT "chr6.pos\tabase\tcbase\tgbase\ttbase\tnbase\tdelte\tinsert\n";

my %depth;

#SAM file 
my $f ="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/5YH.best.sorted.bam";

open IN , "samtools view $f | " and print "$f\n";
my $line = 0;

while (<IN>){
    my %seq_hash;#sequence based (start from 1) hash {ref_pos}{base N.O}{ref}{alt}
    next if (/^@/);
    my ($chr,$pos,$cigar,$seq)=(split /\t/)[2,3,5,9];
    next if ($chr ne "chr6");
    next if ($pos < 28000000);
    last if ($pos > 33500000);

    #$line++;
    #last if ($line>100);
    #print "$line\n";

    $_=~/MD:Z:(.*?)\t/;
    my $md="$1";
    my $now_ci=$cigar;
    my $now_md=$md;
    my $now_pos=$pos;#position in reference genome
    my $now_seq=$seq;
    my $now_base=0;#N.O. of base in sequence

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
                $seq_hash{$now_pos}{$now_base}{$base}{"M"}="mis/match";
            }
        }elsif ($type eq "I"){
            for my $i (1 .. $num){
                $now_base++;
                $now_pos+=0;
                my $base=substr($now_seq,0,1,"");
                $seq_hash{$now_pos}{$now_base}{$base}{"I"}="insert";
            }
        }elsif ($type eq "D"){
            for my $i (1 .. $num){
                $now_base+=0;
                $now_pos++;
                $seq_hash{$now_pos}{$now_base}{"NA"}{"D"}="delet";
            }
        }elsif ($type eq "S"){
            for my $i (1 .. $num){
                $now_base++;
                $now_pos+=0;
                my $base=substr($now_seq,0,1,"");#base of sequenced seq
                $seq_hash{$now_pos}{$now_base}{"NA"}{"S"}="soft";
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
                            for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                                $seq_hash{$krp}{$ksp}{$kb}{$kt}="$kb" if ($kt eq "M");
                                $seq_hash{$krp}{$ksp}{$kb}{$kt}="NA" if ($kt eq "I");
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
                        for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                            if ($kt eq "M"){
                                delete $seq_hash{$krp}{$ksp}{$kb}{"M"};
                                $seq_hash{$krp}{$ksp}{$kb}{"Mis"}="$ref"; 
                            }
                            print"Script error:$ksp $krp $kt == M\n" if($kt eq "I");
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
        for my $ksp (sort {$a<=>$b} keys %{$seq_hash{$krp}}){
            for my $kb (keys %{$seq_hash{$krp}{$ksp}}){
                for my $kt (keys %{$seq_hash{$krp}{$ksp}{$kb}}){
                    my $kr=$seq_hash{$krp}{$ksp}{$kb}{$kt};
                    #print"Script error:$krp\t$ksp\t$kb\t$kt\t$kr\n" if ($kt eq "Mis" && $kb eq $kr);
                    #
                    if($kt eq "D" || $kt eq "I"){
                        $depth{$krp}{$kt}++;
                    }else{
                        $depth{$krp}{$kb}++;
                    }
                }
            }
        }
    }
}
close IN ;

for my $kp (sort {$a<=>$b} keys %depth){
    print "$kp\n";
    my $abase = 0;
    my $tbase = 0;
    my $cbase = 0;
    my $gbase = 0;
    my $nbase = 0;
    my $delet = 0;
    my $insert = 0;
    $abase = $depth{$kp}{"A"} if (exists $depth{$kp}{"A"});
    $cbase = $depth{$kp}{"C"} if (exists $depth{$kp}{"C"});
    $gbase = $depth{$kp}{"G"} if (exists $depth{$kp}{"G"});
    $tbase = $depth{$kp}{"T"} if (exists $depth{$kp}{"T"});
    $nbase = $depth{$kp}{"N"} if (exists $depth{$kp}{"N"});
    $delet = $depth{$kp}{"D"} if (exists $depth{$kp}{"D"});
    #$insert = $depth{$kp}{"I"} if (exists $depth{$kp}{"I"});
    $insert = 1 if (exists $depth{$kp}{"I"});# maybe a lot base insert but only count once
    print OUT "$kp\t$abase\t$cbase\t$gbase\t$tbase\t$nbase\t$delet\t$insert\n";
}
close OUT;
