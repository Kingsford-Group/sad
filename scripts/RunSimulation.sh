#!/bin/bash

codedir=$0
codedir=${codedir%/*}

prepdir=$1

Type=("PC" "Full")
GTFfiles=("${prepdir}/gencode.v26.annotation.pc.gtf" "${prepdir}/gencode.v26.annotation.gtf")
TransFastafiles=("${prepdir}/gencode.v26.pc.transcripts.fa" "${prepdir}/gencode.v26.transcripts.fa")
NumEvents=(200 500 1000 1500)
RefQuantPrefix=("${folder}/../data/BaseExpression/ERR188297" "${folder}/../data/BaseExpression/SRR3192396" "${folder}/../data/BaseExpression/SRR3192412")
RefQuantID=(GEU GM12878 K562)
SADFolder="sad"

#i=0
for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	for n in ${NumEvents[@]}; do
		for ((j=0; j<${#RefQuantPrefix[@]}; j++)); do
			ID=${RefQuantID[${j}]}

			echo ${prepdir}/simu_${t}_${ID}_${n}

			# generate annotation GTF for star alignment
			if [[ ! -e ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf ]]; then
				echo -e "\tWriting new annotation GTF file."
				python3 ${codedir}/ExtractPCAnnotation.py ${prepdir}/simu_${t}_${ID}_${n}_reference.fa ${gtffile} ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf
			fi

			# Aligning reads to reference genome
			if [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam ]]; then
				echo -e "\tAligning reads to reference genome"
				mkdir -p ${prepdir}/Star/${t}_${ID}_${n}/
				STAR --runThreadN 4 --genomeDir ${prepdir}/StarIndex/ --readFilesIn ${prepdir}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${prepdir}/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --readFilesCommand gunzip -c --outFileNamePrefix ${prepdir}/Star/${t}_${ID}_${n}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --sjdbGTFfile ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf --limitBAMsortRAM 32416692217
			fi

			# transcriptome assembly: stringtie
			if [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/stringtiegenes.gtf ]]; then
				echo -e "\tstringtie assembly"
				stringtie ${prepdir}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -G ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf -o ${prepdir}/Star/${t}_${ID}_${n}/stringtiegenes.gtf
			fi

			# transcriptome assembly: scallop
			if [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/scallopgenes.gtf ]]; then
				echo -e "\tscallop assembly"
				scallop --verbose 0 -i ${prepdir}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -o ${prepdir}/Star/${t}_${ID}_${n}/scallopgenes.gtf
			fi

			# summarize assembly result
			if [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms.txt ]] || [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms.txt ]]; then
				echo -e "\tsummarizing transcript assembly result for stringtie"
				gffcompare -o ${prepdir}/Star/${t}_${ID}_${n}/stringtiegffcomp -r ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/Star/${t}_${ID}_${n}/stringtiegenes.gtf
				${codedir}/../bin/assemblypost ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/Star/${t}_${ID}_${n}/stringtiegenes.gtf ${prepdir}/Star/${t}_${ID}_${n}/stringtiegffcomp.stringtiegenes.gtf.tmap ${prepdir}/Star/${t}_${ID}_${n}/stringtie
				echo -e "\tsummarizing transcript assembly result for scallop"
				gffcompare -o ${prepdir}/Star/${t}_${ID}_${n}/scallopgffcomp -r ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/Star/${t}_${ID}_${n}/scallopgenes.gtf
				${codedir}/../bin/assemblypost ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/Star/${t}_${ID}_${n}/scallopgenes.gtf ${prepdir}/Star/${t}_${ID}_${n}/scallopgffcomp.scallopgenes.gtf.tmap ${prepdir}/Star/${t}_${ID}_${n}/scallop
			fi

			# evaluation assembly result
			if [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/evaluation_scallop_multiexon_novelisoforms_existjunc ]] || [[ ! -e ${prepdir}/Star/${t}_${ID}_${n}/evaluation_stringtie_multiexon_novelisoforms_existjunc ]]; then
				echo -e "\tevaluating transcript assembly result for stringtie"
				python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms.txt ${prepdir}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms_uniqgene.txt
				python3 ${codedir}/Simulation.py evatrans ${prepdir}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/evaluation_stringtie_multiexon_novelisoforms_existjunc
				python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms.txt ${prepdir}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms_uniqgene.txt
				python3 ${codedir}/Simulation.py evatrans ${prepdir}/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ${prepdir}/groundtruth_${t}_${ID}_${n}/fusion_noveljunc.txt ${prepdir}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/evaluation_stringtie_noveljunc_novelisoforms_noveljunc

				echo -e "\tevaluating transcript assembly result for scallop"
				python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms.txt ${prepdir}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms_uniqgene.txt
				python3 ${codedir}/Simulation.py evatrans ${prepdir}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/evaluation_scallop_multiexon_novelisoforms_existjunc
				python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms.txt ${prepdir}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms_uniqgene.txt
				python3 ${codedir}/Simulation.py evatrans ${prepdir}/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ${prepdir}/groundtruth_${t}_${ID}_${n}/fusion_noveljunc.txt ${prepdir}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/Star/${t}_${ID}_${n}/evaluation_scallop_noveljunc_novelisoforms_noveljunc
			fi

			# salmon indexing
			if [[ ! -e ${prepdir}/IndexSalmon/${t}_${ID}_${n}/hash.bin ]]; then
				echo -e "\tsalmon indexing"
				mkdir -p ${prepdir}/IndexSalmon/${t}_${ID}_${n}/
				salmon index -t ${prepdir}/simu_${t}_${ID}_${n}_reference.fa -i ${prepdir}/IndexSalmon/${t}_${ID}_${n}/
			fi

			# salmon quantification
			if [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/quant.sf ]]; then
				echo -e "\tsalmon quantification"
				mkdir -p ${prepdir}/salmon/${t}_${ID}_${n}/
				salmon quant -p 4 -l A -i ${prepdir}/IndexSalmon/${t}_${ID}_${n}/ -1 ${prepdir}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz -2 ${prepdir}/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${prepdir}/salmon/${t}_${ID}_${n}/ --writeMappings=${prepdir}/salmon/${t}_${ID}_${n}/mapping.sam
				samtools view -Shb ${prepdir}/salmon/${t}_${ID}_${n}/mapping.sam -o ${prepdir}/salmon/${t}_${ID}_${n}/mapping.bam
				gunzip ${prepdir}/salmon/${t}_${ID}_${n}/aux_info/*gz
			fi

			# expression difference: comparing simulated reads info and salmon quant
			if [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/expression_difference.txt ]]; then
				echo -e "\tcomparing simulated ground truth with salmon quant"
				python3 ${codedir}/GroupTruthExp.py ${prepdir}/simu_${t}_${ID}_${n}_target.fa ${prepdir}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${prepdir}/salmon/${t}_${ID}_${n}/quant.sf ${gtffile} ${prepdir}/salmon/${t}_${ID}_${n}/
			fi

			# run SADpipe
			if [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall ]]; then
				python3 ${codedir}/SADpipe.py -t ${prepdir}/simu_${t}_${ID}_${n}_reference.fa -a ${prepdir}/simu_${t}_${ID}_${n}_annotation.gtf -s ${prepdir}/salmon/${t}_${ID}_${n}/ -o ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad
			fi

			# evaluate SAD new isoform prediction
			if [[ -e ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted ]] && [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_overall_exist ]]; then
				echo -e "\tevaluating SAD new isoform prediction"
				python3 ${codedir}/Pred_Trans2Gene.py ${gtffile} ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted_uniqgene
				python3 ${codedir}/Simulation.py evatrans ${prepdir}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted_uniqgene ${gtffile} ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_overall_exist
			fi

		done
	done
done
