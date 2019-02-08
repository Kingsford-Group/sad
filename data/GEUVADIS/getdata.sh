#!/bin/bash

folder=$0
folder=${folder%/*}

while read -r line; do
	read -ra x <<< "${line}"
	c=${#x[@]}

	if [[ ! -e ${folder}/${x[$((c-1))]}/${x[$((c-1))]}_1.fastq.gz ]]; then
		mkdir -p ${folder}/${x[$((c-1))]}
		wget -P ${folder}/${x[$((c-1))]} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${x[$((c-1))]:0:6}/${x[$((c-1))]}/${x[$((c-1))]}_1.fastq.gz
	fi

	if [[ ! -e ${folder}/${x[$((c-1))]}/${x[$((c-1))]}_2.fastq.gz ]]; then
		mkdir -p ${folder}/${x[$((c-1))]}
		wget -P ${folder}/${x[$((c-1))]} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${x[$((c-1))]:0:6}/${x[$((c-1))]}/${x[$((c-1))]}_2.fastq.gz
	fi
done < ${folder}/Metadata.txt
