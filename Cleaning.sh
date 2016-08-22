#!/bin/bash

# This is a pre-processing script.
# This script takes 4 arguments. 
# First argument is R1_combined.fastq
# Second argument is R2_combined.fastq 
# Third argument is output prefix
# Forth argument is number of threads
# Use following command to run this script
# For example: bash Cleaning.sh R1_combined.fastq R2_combined.fastq output_prefix

usage=`echo -e "\n Usage: Cleaning.sh file_R1.fastq file_R2.fastq OutputPrefix NumThreads \n"`;

if [[ ! $1 ]]
then
        printf "${usage}\n\n";
exit;
fi

echo "Processing Input files" $1 $2

file1=$1
file2=$2

echo $file1;
echo $file2;

if [[ $file1 == *.fastq || $file1 == *.fq ]]
	then
	f1="${file1%.*}";
fi
echo $f1;

if [[ $file2 == *.fastq || $file2 == *.fq ]]
	then
	f2="${file2%.*}";
fi

echo $f2;

#Set location of taxonomy DB
tax_loc=/software/KronaTools-2.6.1/taxonomy/

#set number of threads for the pipeline
threads=$4

#if [[ -z "$4" ]]; then threads = 8; else threads = $4 fi
	

# Count number of reads in given fastq files

touch $3.log
echo "This script started at `date` " > $3.log
echo "Raw number of reads in the fastq files" >> $3.log
echo $1 "contains following number of reads" >> $3.log
expr `(wc -l $1 |cut -f1 -d " ")` / 4 >> $3.log
echo $2 "contains following number of reads" >> $3.log
expr `(wc -l $2 |cut -f1 -d " ")` / 4 >> $3.log

echo "Running trim galore"

#Next step is to remove adapters using Trim Galore

trim_galore --length 80 --paired $1 $2 >/dev/null 2>&1

#Next step is to rename output files from trim_galore

mv ${f1}_val_1.fq $3_val_1.fq
mv ${f2}_val_2.fq $3_val_2.fq

#Count number of reads in given fastq files after running trim_galore

echo "Number of reads after trim_galore" >> $3.log
expr `(wc -l $3_val_1.fq |cut -f1 -d " ")` / 4 >> $3.log
echo "Number of reads after trim_galore" >> $3.log
expr `(wc -l $3_val_2.fq |cut -f1 -d " ")` / 4 >> $3.log

echo "Finish trim Galore now running Prinseq"

#Prinseq - Remove short sequences
 
prinseq-lite.pl -fastq $3_val_1.fq -fastq2 $3_val_2.fq -min_len 50 -out_format 3 -out_good $3_prinseq -out_bad null > /dev/null 2>&1

#Count number of reads in given fastq files after running prinseq

echo $1 "contains following number of reads after prinseq" >> $3.log
expr `(wc -l $3_prinseq_1.fastq|cut -f1 -d " ")` / 4 >> $3.log
echo $2 "contains following number of reads after prinseq" >> $3.log
expr `(wc -l $3_prinseq_2.fastq|cut -f1 -d " ")` / 4 >> $3.log

echo "Running ribopicker to clean against rRNA database"

ribopicker.pl -c 80 -i 90 -f $3_prinseq_1.fastq -dbs slr123,ssr123 -t $threads -id $3_prinseq_1 -keep_tmp_files
ribopicker.pl -c 80 -i 90 -f $3_prinseq_2.fastq -dbs slr123,ssr123 -t $threads -id $3_prinseq_2 -keep_tmp_files

echo "Number of reads with rrna in" $3_prinseq_1.fastq >> $3.log
expr `(wc -l $3_prinseq_1_rrna.fq|cut -f1 -d " ")` / 4 >> $3.log
echo "Number of reads with rrna for" $3_prinseq_1.fastq >> $3.log
expr `(wc -l $3_prinseq_2_rrna.fq|cut -f1 -d " ")` / 4 >> $3.log

echo "Number of reads without rrna in" $3_prinseq_1.fastq >> $3.log
expr `(wc -l $3_prinseq_1_nonrrna.fq|cut -f1 -d " ")` / 4 >> $3.log
echo "Number of reads without rrna for" $3_prinseq_2.fastq >> $3.log
expr `(wc -l $3_prinseq_2_nonrrna.fq|cut -f1 -d " ")` / 4 >> $3.log


mv $3_prinseq_1_nonrrna.fq $3_nonrrna_R1.fq
mv $3_prinseq_2_nonrrna.fq $3_nonrrna_R2.fq

echo "Number of reads submitted to Diamond for " $3_nonrrna_R1.fq >> $3.log
expr `(wc -l $3_nonrrna_R1.fq|cut -f1 -d " ")` / 4 >> $3.log
echo "Number of reads submitted to Diamond for " $3_nonrrna_R2.fq >> $3.log
expr `(wc -l $3_nonrrna_R2.fq|cut -f1 -d " ")` / 4 >> $3.log


#Run Diamond for each fastq files

echo "Running Diamond for" $3_nonrrna_R1.fq "and " $3_nonrrna_R2.fq
mkdir $3_temp_dir
diamond blastx -d $DIAMOND_DB/refseq_protein -p $threads -q $3_nonrrna_R1.fq -a $3_R1_diamond -t $3_temp_dir -k 1
diamond blastx -d $DIAMOND_DB/refseq_protein -p $threads -q $3_nonrrna_R2.fq -a $3_R2_diamond -t $3_temp_dir -k 1

diamond view -a $3_R1_diamond.daa -o $3_R1_diamond.m8 
diamond view -a $3_R2_diamond.daa -o $3_R2_diamond.m8 


ktImportBLAST $3_R1_diamond.m8 -o $3_R1_diamond_krona.html > /dev/null 2>&1
ktImportBLAST $3_R2_diamond.m8 -o $3_R2_diamond_krona.html > /dev/null 2>&1
# 
#Download of NCBI's taxonomy database
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz

#Uncompress gi_taxid_prot.dmp.gz 
#gunzip gi_taxid_prot.dmp.gz
#use dos2unix on the gi_taxid file  just to check formatting

echo "Mapping GI from diamond output to taxonomy database"

#Extraction of GIs from Diamond tabular output in both files

cut -d "|" -f2 $3_R1_diamond.m8 > $3_R1_gi
cut -d "|" -f2 $3_R2_diamond.m8 > $3_R2_gi
echo $3_R1_gi "contains following number of lines after diamond" >> $3.log
wc -l $3_R1_gi >> $3.log
echo $3_R2_gi "contains following number of lines after diamond" >> $3.log
wc -l $3_R2_gi >> $3.log

# Sort and Deduplicate the extract of GIs from both files

sort -n $3_R1_gi | uniq > $3_R1_gi_unique 
sort -n $3_R2_gi | uniq > $3_R2_gi_unique
echo $3_R1_gi_unique "contains following number of unique lines after diamond" >> $3.log
wc -l $3_R1_gi_unique >> $3.log
echo $3_R2_gi_unique "contains following number of unique lines after diamond" >> $3.log
wc -l $3_R2_gi_unique >> $3.log

# Search for patterns (GIs) in both files 
awk 'FNR==NR{a[$1]=1;}FNR!=NR{if(a[$1]==1) print}'  $3_R1_gi_unique $tax_loc/gi_taxid_prot.dmp > $3_R1_match_gi_taxid

awk 'FNR==NR{a[$1]=1;}FNR!=NR{if(a[$1]==1) print}'  $3_R2_gi_unique $tax_loc/gi_taxid_prot.dmp > $3_R2_match_gi_taxid

echo $3_R1_match_gi_taxid "contains following number of GI that has a matching taxa ID" >> $3.log
wc -l $3_R1_match_gi_taxid >> $3.log
echo $3_R2_match_gi_taxid "contains following number of GI that has a matching taxa ID" >> $3.log
wc -l $3_R2_match_gi_taxid >> $3.log

# Extraction tax_id (1st column) and division_id (5th column) from nodes.dmp 
awk -F '\t[|]\t' '{print $1"\t"$5}' $tax_loc/nodes.dmp >$3_nodes_tax_division_id.dmp 

# Cutting and sorting the 2nd column - (only taxid) in both files 
cut -f2 $3_R1_match_gi_taxid > $3_R1_match_taxid
cut -f2 $3_R2_match_gi_taxid > $3_R2_match_taxid
echo $3_R1_match_taxid "contains following number of tax ID" >> $3.log
wc -l $3_R1_match_taxid >> $3.log
echo $3_R2_match_gi_taxid "contains following number of lines of tax ID" >> $3.log
wc -l $3_R2_match_gi_taxid >> $3.log

sort -n $3_R1_match_taxid | uniq > $3_R1_match_taxid_unique
sort -n $3_R2_match_taxid | uniq > $3_R2_match_taxid_unique
echo $3_R1_match_taxid_unique "contains following number of lines of tax ID after sort" >> $3.log
wc -l $3_R1_match_taxid_unique >> $3.log
echo $3_R2_match_taxid_unique "contains following number of lines of tax ID after sort" >> $3.log
wc -l $3_R2_match_taxid_unique >> $3.log

# Search taxid and division in both files 
awk 'FNR==NR{a[$1]=1;}FNR!=NR{if(a[$1]==1) print}' $3_R1_match_taxid_unique $3_nodes_tax_division_id.dmp > $3_R1_taxid_division
awk 'FNR==NR{a[$1]=1;}FNR!=NR{if(a[$1]==1) print}' $3_R2_match_taxid_unique $3_nodes_tax_division_id.dmp > $3_R2_taxid_division

echo $3_R1_taxid_division "contains following number of lines for taxid and division id" >> $3.log
wc -l $3_R1_taxid_division >> $3.log
echo $3_R2_taxid_division "contains following number of lines for taxid and division id" >> $3.log
wc -l $3_R2_taxid_division >> $3.log

echo "Filtering sequences based on the division ID"

# Filtering the rows that has division ids to remove sequences (Bacteria, Invertebrates, Mammals, Phages, Plants, Primates, Rodents, Synthetics, Vertebrates) #both files
awk -F "\t" '{if ($2=="0" || $2=="1" || $2=="2" || $2=="3" || $2=="4"|| $2=="5" || $2=="6" || $2=="7" || $2=="10") print $1"\t"$2}' $3_R1_taxid_division > $3_R1_taxid_to_remove
awk -F "\t" '{if ($2=="0" || $2=="1" || $2=="2" || $2=="3" || $2=="4"|| $2=="5" || $2=="6" || $2=="7" || $2=="10") print $1"\t"$2}' $3_R2_taxid_division > $3_R2_taxid_to_remove

# Sorting and joining files with unique taxaID and file with taxaID and division to get the list of taxaID GI Division information in one file-(output structure is taxid gi division) 

sort -n -k2 $3_R1_match_gi_taxid > $3_R1_match_gi_taxid_sorted
sort -n -k2 $3_R2_match_gi_taxid > $3_R2_match_gi_taxid_sorted

sort -n -k1 $3_R1_taxid_to_remove > $3_R1_taxid_to_remove_sorted
sort -n -k1 $3_R2_taxid_to_remove > $3_R2_taxid_to_remove_sorted

join -1 2 -2 1 $3_R1_match_gi_taxid_sorted $3_R1_taxid_to_remove_sorted > $3_R1_taxid_gi_div
join -1 2 -2 1 $3_R2_match_gi_taxid_sorted $3_R2_taxid_to_remove_sorted > $3_R2_taxid_gi_div

echo $3_R1_taxid_gi_div "contains following number of lines after sorting taxid gi division" >> $3.log
wc -l $3_R1_taxid_gi_div >> $3.log
echo $3_R2_taxid_gi_div "contains following number of lines after sorting taxid gi division" >> $3.log
wc -l $3_R2_taxid_gi_div >> $3.log

# Extract read_id and GIs from diamond output and change "gi" to "nothing" and insert @ in start of name of my reads.
cut -d "|" -f 1,2 $3_R1_diamond.m8| sed -e 's/gi|//g' -e 's/^//' > $3_R1_readid_gi
cut -d "|" -f 1,2 $3_R2_diamond.m8| sed -e 's/gi|//g' -e 's/^//' > $3_R2_readid_gi

echo $3_R1_readid_gi "contains following number of lines after extract read_id and GIs from diamond output" >> $3.log
wc -l $3_R1_readid_gi >> $3.log
echo $3_R2_readid_gi "contains following number of lines after extract read_id and GIs from diamond output" >> $3.log
wc -l $3_R2_readid_gi >> $3.log

# Sorting and joining the readids with the gi, taxid and div
sort -n -k2 $3_R1_readid_gi > $3_R1_readid_gi_sorted
sort -n -k2 $3_R2_readid_gi > $3_R2_readid_gi_sorted

sort -n -k2 $3_R1_taxid_gi_div > $3_R1_taxid_gi_sorted_div
sort -n -k2 $3_R2_taxid_gi_div > $3_R2_taxid_gi_sorted_div

join -1 2 -2 2 $3_R1_readid_gi_sorted $3_R1_taxid_gi_sorted_div > $3_R1_gi_readid_taxid_div
join -1 2 -2 2 $3_R2_readid_gi_sorted $3_R2_taxid_gi_sorted_div > $3_R2_gi_readid_taxid_div

echo $3_R1_gi_readid_taxid_div "contains following number of lines the readids with the gi, taxid and div" >> $3.log
wc -l $3_R1_gi_readid_taxid_div >> $3.log
echo $3_R2_gi_readid_taxid_div "contains following number of lines the readids with the gi, taxid and div" >> $3.log
wc -l $3_R2_gi_readid_taxid_div >> $3.log

# Remove reads from the original cleanup data and extract the column with the readid to a new file needs to be checked
cut -d " " -f2 $3_R1_gi_readid_taxid_div > $3_R1_readid.txt 
cut -d " " -f2 $3_R2_gi_readid_taxid_div > $3_R2_readid.txt 

echo "Removing undesirable sequences"

# Remove sequences from fastq file where list of headers of the sequences are provided in a file - using inverse grep

echo "Filtering fastq files using filter_fastq.pl"

filter_fastq.pl -r $3_R1_readid.txt $3_nonrrna_R1.fq > $3_R1_filtered.fastq 
filter_fastq.pl -r $3_R2_readid.txt $3_nonrrna_R2.fq > $3_R2_filtered.fastq 


#fastqCombinePairedEnd.py $3_R1_filtered.fastq $3_R2_filtered.fastq

prinseq-lite.pl -fastq $3_R1_filtered.fastq -fastq2 $3_R2_filtered.fastq -min_len 10 -out_bad null -out_good $3_clean

#mv $3_R1_filtered.fastq_pairs_R1.fastq $3_clean_R1.fq
#mv $3_R2_filtered.fastq_pairs_R2.fastq $3_clean_R2.fq

echo  "Number of reads in " $3_clean_1.fastq >> $3.log
expr `(wc -l $3_clean_1.fastq|cut -f1 -d " ")` / 4 >> $3.log
echo  "Number of reads in " $3_clean_2.fastq >> $3.log
expr `(wc -l $3_clean_2.fastq|cut -f1 -d " ")` / 4 >> $3.log

rm -fr $3_temp_dir
rm -f $3_R1_taxid* $3_R2_taxid*  $3*.daa $3*unique $3*sorted $3*val* $3*gi $3*txt $3_R1_filtered.fastq $3_R2_filtered.fastq $3_R1_match_taxid $3_R2_match_taxid $3_R1_match_gi_taxid $3_R2_match_gi_taxid $3_prinseq_1_singletons.fastq $3_prinseq_2_singletons.fastq $3_nodes_tax_division_id.dmp

echo "This script finished analysis at `date` " >> $3.log
