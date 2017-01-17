#!/usr/bin/perl
#Run ITSx on all FASTA files in a directory.
#Created 9-1-16 by Bradford Condon
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $usage = "$0 -w directory of blast databases  \n"; 

my $input_file;
my $wd;
#read user options
GetOptions(
	"w=s"  => \$wd
);
if (!$wd){
 $wd = './';
}

mkdir("ITSx");

print $wd;
##list blast databases
my @files = <$wd*.fasta>;
foreach my $file (@files) {   
	my $fileBit = basename($file);
system("ITSx -i $fileBit -t f --cpu 4 -o ITSx/$fileBit.ITS.fasta ");
}




