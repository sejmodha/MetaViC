#!/bin/bash

# This is a De novo assemble script.
# This script takes 3 arguments. 
# First argument is output_clean_R1.fq Second argument is output_clean_R2.fq and Third argument is output prefix
# Use following command to run this script
# For example: bash Assemble.sh output_clean_R1.fq output_clean_R2.fq output_assemble 4

usage=`echo -e "\n Usage: Assemble.sh  file_clean_R1.fastq file_clean_R2.fastq OutputPrefix NumThreads \n"`;

if [[ ! $1 ]]
then
        printf "${usage}\n\n";
exit;
fi

threads=$4

echo "Assemble Input files" $1 $2

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

#De novo assemblies 

touch $3.log

echo "This script started at `date` " > $3.log

echo "Running idba_ud"

fq2fa --merge $1 $2 $3_idba_input.fa 
idba_ud -r $3_idba_input.fa --mink 21 --maxk 121 --step 10 --num_threads $threads -o $3_idba_output 

echo "Running spades"

spades.py -1 $1 -2 $2  -t $threads --only-assembler -o $3_spades_output 

echo "Running GARM"

ln -s $3_idba_output/contig.fa $3_idba_contigs.fa

ln -s $3_spades_output/contigs.fasta $3_spades_contigs.fa

dir=`pwd`

printf $dir'/'$3'_idba_contigs.fa  idba\n'$dir'/'$3'_spades_contigs.fa spades\n' > garm_config

GARM.pl -g garm_config -o $3_contigs_garm

ln -s $3_contigs_garm/final_merged_contigs.fasta $3_garm_contigs.fa

ln -s $3_contigs_garm/final_bin.fasta $3_garm_bin.fa

cat $3_garm_contigs.fa $3_garm_bin.fa > $3_garm_combined.fa

#Align reads back to contigs 

bowtie2-build $3_garm_combined.fa $3_garm_combined

bowtie2 -x $3_garm_combined -1 $1 -2 $2 -S $3.sam -p $threads

samtools view -b $3.sam > $3.bam

samtools sort -o $3_sorted.bam $3.bam 
samtools index $3_sorted.bam

weeSAMv1.1 -b $3_sorted.bam -out $3_weeSAM_out

rm $3.sam $3.bam

bam2fastq --no-aligned --force --strict -o $3_unmapped#.fq $3_sorted.bam
	
#Evaluation of genome assemblies

quast.py -o $3_quast_output  -l "idba_ud, spades, garm" $3_idba_contigs.fa $3_spades_contigs.fa $3_garm_contigs.fa

echo "Running Diamond"

#Run Diamond for GARM contigs

mkdir $3_temp_dir

diamond blastx -d $DIAMOND_DB/refseq_protein -p $threads -q $3_garm_combined.fa -a $3_garm_combined_contigs_diamond -t $3_temp_dir --top 1

diamond view -a $3_garm_combined_contigs_diamond.daa -o $3_garm_combined_contigs_diamond.m8 

#echo "Running PHMMER"

#Run PHMMER for GARM contigs

#phmmer: Search a single protein sequence against a protein sequence database. (BLASTP-like)
#nhmmer: Search a DNA sequence, alignment, or profile HMM against a DNA sequence database. (BLASTN-like)

#phmmer [-options] $3_garm_contigs.fa <seqdb> -o $3_contigs_phmmer

echo "Krona"

#Run Krona

ktImportBLAST $3_garm_combined_contigs_diamond.m8 -o $3_contigs_diamond_krona.html
#ktImportBLAST $3_contigs_phmmer -o $3_contigs_phmmer_krona.html
 
echo "This script finished analysis at `date` " >> $3.log 
 
 
