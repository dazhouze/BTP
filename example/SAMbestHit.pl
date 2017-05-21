#!/usr/bin/perl -w
use strict;
open IN , "$ARGV[0]";
open OUT , ">$ARGV[1]";
my %hash;#store read align
while (<IN>){
    print OUT "$_" if (/^@/);
    next if (/^@/);
    chomp;
    $_=~/AS:i:(.*?)\s/;
    my $as=$1;#AS field
    my ($qname)=(split /\t/)[0];#QNAME to identify read
    $hash{$qname}{$as}="$_";
}
close IN;
#find best hit align
for my $kq (keys %hash){
    my @AS=reverse sort {$a<=>$b} keys %{$hash{$kq}};
       print OUT "$hash{$kq}{$AS[0]}\n";
}
close OUT;
