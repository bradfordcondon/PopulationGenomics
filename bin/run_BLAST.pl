#blasts a single FASTA file against a directory of blast databases.
#output is piped to a single file (.blast extension)
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $usage = "$0 -i Fasta to BLAST as query -w directory of blast databases  \n"; 

my $input_file;
my $wd;
#read user options
GetOptions(
	"i=s"	=> \$input_file,
	"w=s"  => \$wd
);

if (!$wd){
my $wd = './';
}
mkdir("blastReports");
##list blast databases
my @files = <$wd*.fasta>;
foreach my $file (@files) {   
		my $fileBit = basename($file);
system("blastn -query $input_file -out blastReports/$fileBit.blast -db $file -max_target_seqs 1 -outfmt \"6 qseqid length pident sseq\" ");
}
