#blasts a single FASTA file against a directory of blast databases.
#output is piped to a single file (.blast extension)
use strict;
use warnings;
use Getopt::Long;

my $usage = "$0 -i Fasta to BLAST -w directory of blast databases  \n"; 
my $input_file;
my $wd;
#read user options
GetOptions(
	"i=s"	=> \$input_file,
	"w=s"  => \$wd
);
if (!$input_file){
die $usage;
}
if (!$wd){
	print "no blast directory specified, using default\n";
my $wd = '/Users/chet/uky/popGenomeProject7-29-16/wholeGenome/';
}
##list blast databases
my @files = <$wd*.fasta>;
foreach my $file (@files) {   
	my $fileBit = $file;
		$fileBit =~ s/\/Users\/chet\/uky\/popGenomeProject7-29-16\/wholeGenome\///g; #remove path from name file
system("blastn -query $input_file -out blastReports/$fileBit.blast -db $file -max_target_seqs 2 -outfmt \"6 qseqid length pident sseq\" ");
#system("blastn -query $input_file -db $file   >> $input_file.blast");
}




