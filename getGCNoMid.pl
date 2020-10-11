#!/usr/bin/env perl
use strict;
use warnings;
# takes bedtools getfasta "tsv" and computes the %GC
# ignoring the middle base
while (<>){
    chomp;
    my @s = split /\t/;
    my $bases = uc $s[-1];
    my $len = length $bases;
    my $mid = int( $len/2);

    my $numGC = 0;
    my $tot=0;
    for (my $i=0; $i < $len; $i++) {
        if ($i != $mid) {
            my $c = substr($bases, $i, 1);
            if ($c eq 'G' || $c eq 'C') {
                $numGC++;
                $tot++;
            } elsif ($c eq 'A' || $c eq 'T') {
                $tot++;
            }
        }
    }
    my ($chrom, $b, $e) = split /[:-]/, $s[0];
    $mid = int(($e - $b)/2)+1 + $b;
    print $chrom , "\t" , $mid-1 , "\t" , $mid , "\t" , (($e-$b)-1)/2 , "\t" ,  $numGC/$tot , "\n";
       
}


