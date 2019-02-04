#!/bin/bash

folder=$0
folder=${folder%/*}

NR=0

for ID in $(cut -d$'\t' -f10 ${folder}/HumanBodyMap_SraRunTable.txt); do

	((NR++))
	if ((NR==1)); then
		continue
	fi

	if [[ ! -e ${ID}_1.fastq.gz ]]; then
		wget -P ${folder} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ID:0:6}/${ID}/${ID}_1.fastq.gz
	fi

	if [[ ! -e ${ID}_2.fastq.gz ]]; then
		wget -P ${folder} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ID:0:6}/${ID}/${ID}_2.fastq.gz
	fi
done
