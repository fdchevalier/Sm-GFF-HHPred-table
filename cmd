# Date: 2015-05-22

# Download the latest protein dump
dl ftp://ftp.sanger.ac.uk/pub/pathogens/Schistosoma/mansoni/Latest_assembly_annotation_others/Smansoni_chado_dump_protein_20150209.fa

# Split the dump file in individual files
mkdir data
mv Smansoni_chado_dump_protein_20150209.fa data/
cd data
splitfasta.pl Smansoni_chado_dump_protein_20150209.fa
cd ..
mv data/Smansoni_chado_dump_protein_20150209.fa .


# Date: 2015-05-23

# Create list file of sequence to treat and split for parallelization
ls -1 data/* > list 
split -l 100 -a 3 -d list list.d/list.

# Run jobs in parallele
for i in $(ls -1 list.d/*) ; do pqsub hhpred-ann.sh "$i" ; done
