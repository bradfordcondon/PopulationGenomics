#!/usr/bin/env perl
#Random locus generator, version two.
#Takes a BLAST file of reference loci vs all genomes.
#Filters out results not found in all genomes.
#outputs results for each input query as individual FASTA files.

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::searchIO;
use Data::Dumper;

my $directory = "sequences";
 unless(-e $directory or mkdir $directory) {
        die "Unable to create $directory\n";
    }
#global values
my %seqHash;
my %dupeHash;
my %countingHash;
my @loci_list;
my $blastDirectory = '/Users/chet/uky/popGenomeOrthologs9-28-16/supporting/blastReportsv2/';


#read in each blast output file
my @files = <$blastDirectory*>; 

#set 
my $minNumberHits = scalar(@files);

foreach my $file (@files){
	my $prevID = 0;
	my $plength = 0;

open(my $fh, "<", $file)
or die "couldnt open '$file' $!";		
	while (<$fh>){
	chomp;
	my @split = split(/\s+/);#qseqid length pident sseq
	my $qseqid = $split[0];
	my $length = $split[1];
	my $pident = $split[2];
	my $sseq = $split[3];
		if ($length > $plength){
		$seqHash{$file}{$qseqid}{'sequence'} = $sseq;
		$seqHash{$file}{$qseqid}{'length'} = $length;					
		$seqHash{$file}{$qseqid}{'pident'} = $pident;	
		$countingHash{$qseqid}++;								
		}
	$plength = $length;
	}
	
close $fh;
}


##screen out sequences with more than 1 hit above a certain threshhold

#looks like 76 is hte reasonable number here.

	###Build a list of loci to shuffle, but leave out those not found in at least 1 genome except DSLiz!!
	foreach my $k (keys %countingHash) {  
            if (exists $dupeHash{$k}) {  
       		}
                else {
                  	#check for number of keys
                  		if ($countingHash{$k} >= $minNumberHits) {
							push (@loci_list, $k);
						}
				}
		}


	foreach my $locus (@loci_list){
		my $namefile = $locus;
		 $namefile =~ s/scaffold/s/g;

		my $output_file = "MPG_all.fasta";

	open( OUT, ">$output_file" )
	or die "error: cannot open $output_file for writing: $!\n";

		foreach my $database (sort keys %seqHash) {  
									
						my $seqfinal = $seqHash{$database}{$locus}{"sequence"};
						print  OUT ">$database\n$seqfinal\n";
		}
			
	}

close OUT;

