#
#Bradford Condon
#9-30-16
#Goal of script:

#*read in ITS1 and ITS2 FASTAs seperately
#* Exclude regions less than (75?) basepairs, greater than (400?) basepairs.  (Some outputs have very very large regions included- looks like the entire scaffold, actually)
#Output two FASTA files, one for ITS1 and one for ITS2

use strict;
use warnings;
use Bio::SeqIO;


my $itsxDirectory = './';
my %ITS1hash;
my %ITS2hash;


my @files = <$itsxDirectory*.fasta>;
foreach my $file (@files) {  

	if ($file  =~ /.ITS1.fasta/){ #matches ITS1 

		my $input = Bio::SeqIO->new (-file => "<$file", '-format' => 'Fasta')
			or die "error: failure opening fasta $file for reading: $!\n";

		while( my $seqObject = $input->next_seq ){
			my $id  = $seqObject->id;
			my $seq = $seqObject->seq;
			my $length = $seqObject->length;
			if ($length >= 70 && $length <= 400){
				$ITS1hash{$id} = $seq;
			}
		}
	}

	if ($file  =~ /.ITS2.fasta/){
		
		my $input = Bio::SeqIO->new (-file => "<$file", '-format' => 'Fasta')
			or die "error: failure opening fasta $file for reading: $!\n";

		while( my $seqObject = $input->next_seq ){
			my $id  = $seqObject->id;
			my $seq = $seqObject->seq;
			my $length = $seqObject->length;
			if ($length >= 70 && $length <= 400){
				$ITS2hash{$id} = $seq;
			}
		}
	}
}

my $output_file = "ITS1extracted.fasta";
open( OUT, ">$output_file" )
	or die "error: cannot open $output_file for writing: $!\n";


foreach my $k (keys %ITS1hash) {
	my $outSeq = $ITS1hash{$k};
	print OUT ">","$k", "\n", $outSeq, "\n";
}
close OUT;


my $output_file2 = "ITS2extracted.fasta";
open( OUT2, ">$output_file2" )
or die "error: cannot open $output_file2 for writing: $!\n";

foreach my $k (keys %ITS2hash) {
	my $outSeq = $ITS2hash{$k};
	print OUT2 ">","$k", "\n", $outSeq, "\n";
}
close OUT2;


