#!/bin/bash
#SBATCH -J longreadASE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=56000M
#SBATCH --time=5:00:00
#SBATCH --export=none
#SBATCH --output=./slurm-%j.out
##****************************##
##    ASE longread workflow   ##
##****************************##
## Project: longreads
## Lison Lemoine
## Creation: 10.2023
## Submitting bam files to separate bam script

set -e

error_message="
!! An environment must be created for whatshap !! 

Use the following command:                        
conda create -n whatshap-17 whatshap=1.7          

Requirements:
	samtools=1.12
	pbmm2=1.10
	htseq=0.12.4
"
date
echo "longreads ASE"

module load python/2.7_intel
source activate whatshap-17 || echo $error_message

######------ Set environment ------######
while getopts ":c:" flag; do
	case $flag in
	c) configfile=$OPTARG;;
esac
done

#include config file
if [ ! -s $configfile ]; then
	echo "Config file not found or no annotation file given."
	exit 4
fi

params=("pipe_location" "ref" "name" "vcf" "annotation" "aligned")

for param in "${params[@]}" 
do
	param_line=`grep -v "\#" ${configfile} | grep "$param"`
	eval "$param_line"
done

input=${pipe_location}"/input"
OUTdir=${pipe_location}"/output"

######------ Checkers ------######
# check if flag parameters are set
if [ ! -s $OUTdir ]; then
	echo "
Error: output folder ${OUTdir} can not be found."
	exit 1
fi
if [ ! -s $ref ]; then
	echo "
Error: reference file ${ref} can not be found."
	exit 2
fi
if [ ! -s $vcf ]; then
	echo "
Error: vcf file ${vcf} can not be found.
Please generate file using the create_vcf.sh script in the supplemental_scripts folder"
	exit 3
fi

echo "pipeline location: $pipe_location"
echo "outputdir: $OUTdir"
echo "annotation: $annotation"
echo "reference: $ref"
echo "snp_file: $vcf"

######------ Workflow ------######
# the following steps are done via SMRTLink so they can be ignored
# also the barcodes and primer fasta file needs to be known
# lima hifi_reads.bam barcoded_primers.fasta fl.bam --isoseq --peek-guess
# isoseq refine fl.bam primers.fasta flnc.bam

## ALIGN
if [ -z $aligned ]; then
	aligned=${input}"/flnc.aligned.bam"

	if [ ! -s $aligned ]; then
	echo "Did not found an aligned file, generating one"
    bam=$input/${name}_flnc.bam

   		if [ ! -s $bam ]; then
		echo "
		Error: Input BAM file ${bam} can not be found."
		exit 3
    	fi

    pbmm2 align $ref $bam $aligned --sort -j 4 -J 2 --log-level INFO
	fi
fi

## TAG isoforms
echo "tagging isoforms using whatshap"

whatshap haplotag -o ${OUTdir}/${name}_haplotagged.bam --reference $ref $vcf $aligned --ignore-read-groups --ignore-linked-read --skip-missing-contigs --output-haplotag-list $OUTdir/${name}_haplotypes.tsv

# to visualize in IGV
samtools index ${OUTdir}/${name}_haplotagged.bam

## split
# this uses the original aligned bam file
# h1 is maternal and h2 is paternal
whatshap split --output-h1 ${OUTdir}/${name}_h1.bam --output-h2 ${OUTdir}/${name}_h2.bam $aligned $OUTdir/${name}_haplotypes.tsv
# to visualize in IGV
samtools index -M ${OUTdir}/${name}_h2.bam ${OUTdir}/${name}_h1.bam

## Counting

samtools sort -m 6G ${OUTdir}/${name}_h1.bam -o $OUTdir/sorted_${name}_h1.bam
samtools sort -m 6G ${OUTdir}/${name}_h2.bam -o ${OUTdir}/sorted_${name}_h2.bam
htseq-count -s yes -r pos -f bam $OUTdir/sorted_${name}_h1.bam ${annotation} > ${OUTdir}/${name}"_h1.count"
htseq-count -s yes -r pos -f bam $OUTdir/sorted_${name}_h2.bam ${annotation} > ${OUTdir}/${name}"_h2.count"
