#!/bin/bash

codedir=$0
codedir=${codedir%/*}

prepdir=$1

MetaFile="${codedir}/../data/HumanBodyMap/HumanBodyMap_SraRunTable.txt"

Type=("Full")
GTFfiles=("${prepdir}/gencode.v26.annotation.gtf")
TransFastas=("${prepdir}/gencode.v26.transcripts.fa")
SalmonIndex=("${prepdir}/gencode.v26.full")
ReadFolder="${codedir}/../data/HumanBodyMap"
OutDirectory="${prepdir}/HumanBodyMap"
SADFolder="sad"

count=0

i=0
#for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastas[${i}]}
	salmonindex=${SalmonIndex[${i}]}
	for ID in $(cut -d$'\t' -f10 ${MetaFile}); do
		((count++))
		if ((count<=1)); then
			continue
		fi

		read1=${ReadFolder}/${ID}"_1.fastq.gz"
		read2=${ReadFolder}/${ID}"_2.fastq.gz"

		echo "${OutDirectory}/salmon_${t}_${ID}/"

		# running salmon quantification
		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/quant.sf ]]; then
			mkdir -p ${OutDirectory}/salmon_${t}_${ID}/
			echo -e "\tsalmon quantification"

			salmon quant -p 4 -l A -i ${salmonindex} -1 ${read1} -2 ${read2} --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${OutDirectory}/salmon_${t}_${ID}/ --writeMappings=${OutDirectory}/salmon_${t}_${ID}/mapping.sam
			samtools view -Shb ${OutDirectory}/salmon_${t}_${ID}/mapping.sam -o ${OutDirectory}/salmon_${t}_${ID}/mapping.bam
		fi

		# run SADpipe
		if [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue.tsv ]]; then
			python3 ${codedir}/SADpipe.py -t ${transfasta} -a ${gtffile} -s ${OutDirectory}/salmon_${t}_${ID}/ -o ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad
		fi

		# post-process SAD
		if [[ -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue.tsv ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ]]; then
			echo -e "\tpost-process SAD p values"
			python3 ${codedir}/SortPValue.py 1 ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv
			python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv
		fi

		# STAR alignment for transcriptome assembly
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam ]]; then
			mkdir -p ${OutDirectory}/star_${t}_${ID}/
			echo -e "\tSTAR aligning"
			STAR --runThreadN 4 --genomeDir ${prepdir}/StarIndex/ --readFilesIn ${read1} ${read2} --readFilesCommand gunzip -c --outFileNamePrefix ${OutDirectory}/star_${t}_${ID}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --sjdbGTFfile ${gtffile}
		fi

		# transcriptome assembly: Scallop
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ]]; then
			echo -e "\tScallop"
			scallop -i ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam --verbose 0 -o ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf
			gffcompare -o ${OutDirectory}/star_${t}_${ID}/gffcompscallop -r ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf
		fi

		# transcriptome assembly StringTie
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ]]; then
			echo -e "\tStringTie"
			stringtie ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam -G ${gtffile} -o ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf
			gffcompare -o ${OutDirectory}/star_${t}_${ID}/gffcompstringtie -r ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf
		fi

		# analysis: compare SAD and scallop
		if [[ -e ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_pvalue_overall_scallopcomp2 ]]; then
			echo -e "\tcomparing SAD and Scallop prediction"
			python3 ${codedir}/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_pvalue_overall_scallopcomp
		fi

		# analysis: compare SAD and stringtie
		if [[ -e ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_pvalue_overall_stringtiecomp2 ]]; then
			echo -e "\tcomparing SAD and StringTie prediction"
			python3 ${codedir}/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_pvalue_overall_stringtiecomp
		fi

	done
#done


# extracting number of reads per sample
#ReadSumFile=${OutDirectory}/"ReadSummary.txt"
#echo -e "# SampleID\tNumReads" > ${ReadSumFile}

#while read -r line; do
#        t=${Type[${i}]}
#        ID=${line}

#        read -ra y <<< $(grep "total reads" ${OutDirectory}/salmon_${t}_${ID}/logs/salmon_quant.log)
#        nreads=${y[5]}
#        echo -e ${ID}"\t"${nreads} >> ${ReadSumFile}
#done < ${MetaFile}
