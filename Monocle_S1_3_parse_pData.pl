#!/usr/bin/perl

use strict;
use warnings;

my $index = 1;
my %index_RawBatch;
open IN,"<","/Users/m/tmp/data/_AGG11_170901/AGG11_mapped/outs/aggregation_csv.csv" or die $!;
<IN>;
while(<IN>){
	# print "$index\t$_";
	chomp;
	my @in = split/\,/;
	# print "$index\t$in[0]\n";
	$index_RawBatch{$index} = $in[0];
	$index++;
}
close IN;

my %barcode_MEFXENES;
open IN,"<","MEFXENES.csv" or die $!;
<IN>;
while(<IN>){
	# print "$_";
	chomp;
	my @in = split/\,/;
	$barcode_MEFXENES{$in[0]} = $in[1];
}
close IN;

open IN,"<","pData_raw.csv" or die $!;
open OUT,">","pData_to_monocle.csv" or die $!;
<IN>;
print OUT ",Barcode,RawBatch,Batch\n";
while(<IN>){
	# print "$_";
	chomp;
	my @in = split/\,/;
	my $barcode = $in[1];
	# print "$barcode\n";
	my $batchIndex = (split/\-/,$barcode)[1];
	# print "$barcode\t$batchIndex\n";
	if($batchIndex == 1){
		# print "$_\n";
		print OUT "$in[0]\,$barcode\,$index_RawBatch{$batchIndex}\,$barcode_MEFXENES{$barcode}\n";
	}else{
		print OUT "$in[0]\,$barcode\,$index_RawBatch{$batchIndex}\,$index_RawBatch{$batchIndex}\n";
	}
}
close IN;
close OUT;
