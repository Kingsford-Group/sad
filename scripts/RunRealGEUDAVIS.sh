#!/bin/bash

MetaFile="/home/congm1/savanna/savannacong33/RawData/GEUVADIS/Metadata.txt"

Type=("PC" "Full")
GTFfiles=("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.pc.gtf" "/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
SalmonIndex=("/mnt/disk33/user/congm1/NCBI/gencode.v26.pcERCC" "/mnt/disk33/user/congm1/NCBI/gencode.v26.full")
TransFastafiles=("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.pc.transcripts.fa" "/home/congm1/savanna/savannacong33/NCBI/gencode.v26.full.transcripts.fa")
ReadFolder="/home/congm1/savanna/savannacong33/RawData/GEUVADIS"

OutDirectory="/home/congm1/savanna/savannacong33/SADrealdata/GEUVADIS"
SADFolder="test_correctapprox8"

count=0

i=1
#for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastafiles[${i}]}
	salmonindex=${SalmonIndex[${i}]}
	while read -r line; do

		((count++))
		#if ((count<=0)) || ((count>6)); then
		#	continue
		#fi

		#if ((count != 24)) && ((count != 30)); then
		#	continue
		#fi

		read -ra x <<< ${line}
		ID=${x[${#x[@]}-1]}
		read1=${ReadFolder}/${ID}/${ID}"_1.fastq.gz"
		read2=${ReadFolder}/${ID}/${ID}"_2.fastq.gz"

		echo "${OutDirectory}/salmon_${t}_${ID}/"

		# running salmon quantification
		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/quant.sf ]]; then
			mkdir -p ${OutDirectory}/salmon_${t}_${ID}/
			echo -e "\tsalmon quantification"

			salmon quant -p 4 -l A -i ${salmonindex} -1 ${read1} -2 ${read2} --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${OutDirectory}/salmon_${t}_${ID}/ --writeMappings=${OutDirectory}/salmon_${t}_${ID}/mapping.sam
			samtools view -Shb ${OutDirectory}/salmon_${t}_${ID}/mapping.sam -o ${OutDirectory}/salmon_${t}_${ID}/mapping.bam
			gunzip ${OutDirectory}/salmon_${t}_${ID}/aux_info/*gz
		fi

		# expected distribution with bias correction
		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/correction.dat ]]; then
			echo -e "\texpected distribution calculation"
			/home/congm1/savanna/savannacong33/Code/SAD/bin/readsalmonbias correction ${OutDirectory}/salmon_${t}_${ID}/aux_info/ ${transfasta} ${OutDirectory}/salmon_${t}_${ID}/quant.sf ${OutDirectory}/salmon_${t}_${ID}/correction.dat
		fi

		# observed distribution
		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/startpos.dat ]]; then
			echo -e "\tobserved distribution calculation"
			/home/congm1/savanna/savannacong33/Code/SAD/bin/transcovdist2 0 ${OutDirectory}/salmon_${t}_${ID}/quant.sf ${OutDirectory}/salmon_${t}_${ID}/aux_info/eq_classes.txt ${OutDirectory}/salmon_${t}_${ID}/mapping.bam ${OutDirectory}/salmon_${t}_${ID}/startpos.dat
		fi

		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall ]]; then
			echo -e "\tRunning SAD"
			mkdir -p ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/
			/home/congm1/savanna/savannacong33/Code/SAD/bin/SAD ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/quant.sf ${OutDirectory}/salmon_${t}_${ID}/correction.dat ${OutDirectory}/salmon_${t}_${ID}/startpos.dat ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test
		fi

		# post-process SAD
		if [[ -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene ]]; then
			echo -e "\tpost-process SAD p values"
			python3 /home/congm1/savanna/savannacong33/Code/SAD/src/AdjustPValue.py 1 ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted
			python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene
		fi

		## adding DeletionRegion information
		#if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene_append ]]; then
		#	/home/congm1/savanna/savannacong33/Code/SAD/bin/postSAD ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/quant.sf ${OutDirectory}/salmon_${t}_${ID}/correction.dat ${OutDirectory}/salmon_${t}_${ID}/startpos.dat ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test
		#fi

		## Analysis: length rank
		#if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/analysis_lengthrank.txt ]]; then
		#	echo -e "\tAnalysis: length rank"
		#	python3 /mnt/disk33/user/congm1/Code/SAD/src/SADAnalysis_singlesample.py lengthrank ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene_append ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/analysis_lengthrank.txt
		#fi

		# STAR alignment for transcriptome assembly
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam ]]; then
			mkdir -p ${OutDirectory}/star_${t}_${ID}/
			echo -e "\tSTAR aligning"
			STAR --runThreadN 4 --genomeDir /home/congm1/savanna/savannacong33/NCBI/StarIndex_full/ --readFilesIn ${read1} ${read2} --readFilesCommand gunzip -c --outFileNamePrefix ${OutDirectory}/star_${t}_${ID}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15
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
		if [[ -e ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_scallopcomp2 ]]; then
			echo -e "\tcomparing SAD and Scallop prediction"
			#python3 src/FindSADpredinAssembly.py 0 ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene_append ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_scallopcomp
			python3 src/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_scallopcomp2
		fi

		# analysis: compare SAD and stringtie
		if [[ -e ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_stringtiecomp2 ]]; then
			echo -e "\tcomparing SAD and StringTie prediction"
			#python3 src/FindSADpredinAssembly.py 0 ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene_append ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_stringtiecomp
			python3 src/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_sorted_uniqgene ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/test_pvalue_overall_stringtiecomp2
		fi

	done < ${MetaFile}
#done

# extracting number of reads per sample
#ReadSumFile=${OutDirectory}/"ReadSummary.txt"
#echo -e "# SampleID\tNumReads" > ${ReadSumFile}

#while read -r line; do
#	t=${Type[${i}]}
#	read -ra x <<< ${line}
#	ID=${x[${#x[@]}-1]}

#	read -ra y <<< $(grep "total reads" ${OutDirectory}/salmon_${t}_${ID}/logs/salmon_quant.log)
#	nreads=${y[5]}
#	echo -e ${ID}"\t"${nreads} >> ${ReadSumFile}
#done < ${MetaFile}
