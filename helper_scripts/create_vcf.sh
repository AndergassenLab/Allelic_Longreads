#!/bin/bash
##*****************************##
## 	          VCF     	       ##
##*****************************##
## Project: long reads
## Lison Lemoine
## Last modification 10.2023
## Creation: 04.2023
## creates vcf file

date
echo "create vcf"

# requirements
# htslib=1.18

while getopts :s:p: flag; do
    case $flag in
        s) snp_file=$OPTARG;;
        p) pipeline=$OPTARG;;
    esac
done
# check if required parameters are empty -> if so invoke help
[ -z $snp_file ] && exit 1
[ -z $pipeline ] && exit 1

# check if flag parameters are set
if [ ! -s $snp_file ]; then
	echo "
Error: bam folder ${snp_file} can not be found."
	exit 0
fi
if [ ! -s $pipeline ]; then
	echo "
Error: path to workflow ${pipeline} can not be found."
	exit 1
fi

OUTdir=${pipeline}"/input"

#format:
####fileformat=VCFv4.1
##CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
# 1     12  .   T   A  50.0 .      .    GT     0/1 

echo $snp_file

touch ${OUTdir}/snps.vcf
awk '{ print $1 }' $snp_file > ${OUTdir}/1_temp.txt

awk '{ print $2 }' $snp_file > ${OUTdir}/2_temp.txt
#awk -F, '{$1=$1-1}1' ${OUTdir}/2_temp.txt >${OUTdir}/temp && mv ${OUTdir}/temp ${OUTdir}/2_temp.txt
awk -F, '{$1=$1}1' ${OUTdir}/2_temp.txt >${OUTdir}/temp && mv ${OUTdir}/temp ${OUTdir}/2_temp.txt

awk '$4>0 { print substr($4,1,1)}' $snp_file > ${OUTdir}/3_temp.txt
awk '$4>0 { print substr($4,2,2)}' $snp_file > ${OUTdir}/4_temp.txt

paste ${OUTdir}/1_temp.txt ${OUTdir}/2_temp.txt > ${OUTdir}/snps_temp.vcf
awk -F '\t' -v OFS='\t' '{ $(NF+1) = "."; print }' ${OUTdir}/snps_temp.vcf > ${OUTdir}/snps_temp0.vcf
paste ${OUTdir}/snps_temp0.vcf ${OUTdir}/3_temp.txt ${OUTdir}/4_temp.txt > ${OUTdir}/snps_temp.vcf
awk -F '\t' -v OFS='\t' '{ $(NF+1) = "50.0\t.\t.\tGT\t0|1"; print }' ${OUTdir}/snps_temp.vcf > ${OUTdir}/snps_temp0.vcf

( echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"; cat ${OUTdir}/snps_temp0.vcf ) > ${OUTdir}/snps_temp.vcf
( echo -e "##fileformat=VCFv4.1"; cat ${OUTdir}/snps_temp.vcf ) > ${OUTdir}/snps.vcf

rm *temp*

bgzip ${OUTdir}/snps.vcf
tabix ${OUTdir}/snps.vcf.gz
