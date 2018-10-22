#!/usr/bin/perl

use strict;
use warnings;

open IN,"<","fData_raw.csv" or die $!;
open OUT,">","fData_to_monocle.csv" or die $!;
<IN>;
print OUT ",id,gene_short_name\n";
while(<IN>){
	chomp;
	print OUT "$_\n";
}
close IN;
close OUT;
