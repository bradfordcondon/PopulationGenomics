##clean FASTA file names

for f in *.fasta 

do

perl -i -pe 's/\s.*\n/\n/g' $f

done