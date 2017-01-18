#!/usr/bin/env perl
#combine ITS BLASTs
#Created 9-25-16
#last modified 1-17-17
#Takes a BLAST files of ITS vs all genomes.
#Combines into one FASTA file

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::searchIO;
use Data::Dumper;
use File::Basename;



my $directory = "sequences";
#Declare global values
my %seqHash;
my %dupeHash;
my %countingHash;
my @loci_list;
my $blastDirectory;

GetOptions(
	"w=s"  => \$blastDirectory,
	"o=s" => \$namefile
);

if (!$blastDirectory){ #default to current directory
my $blastDirectory = './';
}
if (!$namefile){ #name default output FASTA
	my $namefile = "all_Loci.fasta";
}



#read in each blast output file
my @files = <$blastDirectory*>; 
foreach my $file (@files){
open(my $fh, "<", $file)
or die "couldnt open '$file' $!";		
	while (<$fh>){
	chomp;
	my @split = split(/\s+/);#qseqid length pident sseq
	my $qseqid = $split[0];
	my $length = $split[1];
	my $pident = $split[2];
	my $sseq = $split[3];
	my $fileBit = basename($file);

	$seqHash{$fileBit}{'sequence'} = $sseq;
	}
}

	open( OUT, ">$namefile" )
	or die "error: cannot open $namefile for writing: $!\n";

#print output
	
		foreach my $database (sort keys %seqHash) { 
		print OUT "\n>$database\n"; 
						my $seqfinal = $seqHash{$database}{'sequence'};
						print  OUT "$seqfinal";
	}
close OUT;
