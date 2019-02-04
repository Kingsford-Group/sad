#!/bin/bash

Type=("PC" "Full")
GTFfiles=("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.pc.gtf" "/home/congm1/savanna/savannacong33/NCBI/gencode.v26.annotation.gtf")
TransFastafiles=("/home/congm1/savanna/savannacong33/NCBI/gencode.v26.pc.transcripts.fa" "/home/congm1/savanna/savannacong33/NCBI/gencode.v26.transcripts.fa")
NumEvents=(200 500 1000 1500)
#NumEvents=(1000 1500)
RefQuantFolders=("/home/congm1/savanna/savannacong33/StarAlign/GEUVADIS/ERR188297/" "/home/congm1/savanna/savannacong33/StarAlign/ENCODE/SRR3192396/" "/home/congm1/savanna/savannacong33/StarAlign/ENCODE/SRR3192412/")
RefQuantID=(GEU GM12878 K562)

OutDirectory="/home/congm1/savanna/savannacong33/SADstandardsimu"
SADFolder="test_correctapprox8"

#i=0
for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastafiles[${i}]}
	for n in ${NumEvents[@]}; do
		for ((j=0; j<${#RefQuantFolders[@]}; j++)); do
			ID=${RefQuantID[${j}]}
			if ((i==0)); then
				refquantfile=${RefQuantFolders[${j}]}"salmon/quant.sf"
			else
				refquantfile=${RefQuantFolders[${j}]}"salmon_full/quant.sf"
			fi

			echo ${OutDirectory}/simu_${t}_${ID}_${n}

			# simulating deletion and fusion
			if [[ ! -e ${OutDirectory}/simu_${t}_${ID}_${n}_theoexp.txt ]]; then
				echo -e "\tSimulating deletion and fusion."
				mkdir -p ${OutDirectory}/simu_${t}_${ID}_${n}
				python3 /mnt/disk33/user/congm1/Code/SAD/src/Simulation.py ${transfasta} ${gtffile} ${refquantfile} ${OutDirectory}/simu_${t}_${ID}_${n} ${n} ${n}
			fi

			if [[ ! -e ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ]]; then
				echo -e "\tCategorizing simulated events."
				mkdir -p ${OutDirectory}/groundtruth_${t}_${ID}_${n}
				/home/congm1/savanna/savannacong33/Code/SAD/bin/categorizesimulation ${gtffile} ${OutDirectory}/simu_${t}_${ID}_${n}_removelist.txt ${OutDirectory}/simu_${t}_${ID}_${n}_target.fa ${OutDirectory}/groundtruth_${t}_${ID}_${n}/
			fi

			# simulating reads with polyester
			if [[ ! -e ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ]] && [[ ! -e  ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_1.fasta ]]; then
				echo -e "\tSimulating reads with polyester."
				mkdir -p ${OutDirectory}/polyester_${t}_${ID}_${n}
				/opt/local/stow/R-3.3.2/bin/Rscript /mnt/disk33/user/congm1/Code/SAD/src/SimulationReads.R ${OutDirectory}/simu_${t}_${ID}_${n}_target.fa ${OutDirectory}/simu_${t}_${ID}_${n}_theoexp.txt ${OutDirectory}/polyester_${t}_${ID}_${n}/ 1
			fi

			# gzipping simulated reads
			if [[ $(find ${OutDirectory}/polyester_${t}_${ID}_${n}/ -name *fasta) != "" ]]; then
				echo -e "\tGzipping simulated reads."
				gzip ${OutDirectory}/polyester_${t}_${ID}_${n}/*fasta
			fi

			# generate annotation GTF for star alignment
			if [[ ! -e ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ]]; then
				echo -e "\tWriting new annotation GTF file."
				python3 src/tmpGtfForStar.py ${OutDirectory}/simu_${t}_${ID}_${n}_removelist.txt ${gtffile} ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf
			fi

			# Aligning reads to reference genome
			if [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam ]]; then
				echo -e "\tAligning reads to reference genome"
				mkdir -p ${OutDirectory}/Star/${t}_${ID}_${n}/
				STAR --runThreadN 4 --genomeDir /home/congm1/savanna/savannacong33/NCBI/StarRawIndex/ --readFilesIn ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --readFilesCommand gunzip -c --outFileNamePrefix ${OutDirectory}/Star/${t}_${ID}_${n}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --sjdbGTFfile ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf --limitBAMsortRAM 32416692217
			fi

			# transcriptome assembly: stringtie
			if [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegenes.gtf ]]; then
				echo -e "\tstringtie assembly"
				stringtie ${OutDirectory}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -G ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf -o ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegenes.gtf
			fi

			# transcriptome assembly: scallop
			if [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgenes.gtf ]]; then
				echo -e "\tscallop assembly"
				scallop --verbose 0 -i ${OutDirectory}/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -o ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgenes.gtf
			fi

			# summarize assembly result
			if [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms.txt ]] || [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms.txt ]]; then
				echo -e "\tsummarizing transcript assembly result for stringtie"
				gffcompare -o ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegffcomp -r ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegenes.gtf
				/home/congm1/savanna/savannacong33/Code/SAD/bin/assemblypost ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegenes.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/stringtiegffcomp.stringtiegenes.gtf.tmap ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie
				echo -e "\tsummarizing transcript assembly result for scallop"
				gffcompare -o ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgffcomp -r ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgenes.gtf
				/home/congm1/savanna/savannacong33/Code/SAD/bin/assemblypost ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgenes.gtf ${OutDirectory}/Star/${t}_${ID}_${n}/scallopgffcomp.scallopgenes.gtf.tmap ${OutDirectory}/Star/${t}_${ID}_${n}/scallop
			fi

			# evaluation assembly result
			if [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_scallop_multiexon_novelisoforms_existjunc ]] || [[ ! -e ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_stringtie_multiexon_novelisoforms_existjunc ]]; then
				echo -e "\tevaluating transcript assembly result for stringtie"
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms.txt ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms_uniqgene.txt
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Simulation.py evatrans ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${OutDirectory}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms_uniqgene.txt ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_stringtie_multiexon_novelisoforms_existjunc
				head -200 ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_stringtie_multiexon_novelisoforms_existjunc | cut -f4 | sort | uniq -c
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms.txt ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms_uniqgene.txt
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Simulation.py evatrans ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ${OutDirectory}/groundtruth_${t}_${ID}_${n}/fusion_noveljunc.txt ${OutDirectory}/Star/${t}_${ID}_${n}/stringtie_noveljunc_novelisoforms_uniqgene.txt ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_stringtie_noveljunc_novelisoforms_noveljunc
				head -200 ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_stringtie_noveljunc_novelisoforms_noveljunc | cut -f4 | sort | uniq -c

				echo -e "\tevaluating transcript assembly result for scallop"
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms.txt ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms_uniqgene.txt
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Simulation.py evatrans ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${OutDirectory}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms_uniqgene.txt ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_scallop_multiexon_novelisoforms_existjunc
				head -200 ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_scallop_multiexon_novelisoforms_existjunc | cut -f4 | sort | uniq -c
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms.txt ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms_uniqgene.txt
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Simulation.py evatrans ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ${OutDirectory}/groundtruth_${t}_${ID}_${n}/fusion_noveljunc.txt ${OutDirectory}/Star/${t}_${ID}_${n}/scallop_noveljunc_novelisoforms_uniqgene.txt ${gtffile} ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_scallop_noveljunc_novelisoforms_noveljunc
				head -200 ${OutDirectory}/Star/${t}_${ID}_${n}/evaluation_scallop_noveljunc_novelisoforms_noveljunc | cut -f4 | sort | uniq -c
			fi

			# salmon indexing
			if [[ ! -e ${OutDirectory}/IndexSalmon/${t}_${ID}_${n}/hash.bin ]]; then
				echo -e "\tsalmon indexing"
				mkdir -p ${OutDirectory}/IndexSalmon/${t}_${ID}_${n}/
				salmon index -t ${OutDirectory}/simu_${t}_${ID}_${n}_reference.fa -i ${OutDirectory}/IndexSalmon/${t}_${ID}_${n}/
			fi

			# salmon quantification
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/quant.sf ]]; then
				echo -e "\tsalmon quantification"
				mkdir -p ${OutDirectory}/salmon/${t}_${ID}_${n}/
				salmon quant -p 4 -l A -i ${OutDirectory}/IndexSalmon/${t}_${ID}_${n}/ -1 ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz -2 ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${OutDirectory}/salmon/${t}_${ID}_${n}/ --writeMappings=${OutDirectory}/salmon/${t}_${ID}_${n}/mapping.sam
				samtools view -Shb ${OutDirectory}/salmon/${t}_${ID}_${n}/mapping.sam -o ${OutDirectory}/salmon/${t}_${ID}_${n}/mapping.bam
				gunzip ${OutDirectory}/salmon/${t}_${ID}_${n}/aux_info/*gz
			fi

			# salmon bootstrapping
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/aux_info/bootstrap/quant_bootstraps.tsv ]]; then
				python ~/ocean/oceancong02/Software/salmon-0.8.2/scripts/ConvertBootstrapsToTSV.py ${OutDirectory}/salmon/${t}_${ID}_${n}/ ${OutDirectory}/salmon/${t}_${ID}_${n}/aux_info/bootstrap/
			fi

			# expected distribution with bias correction
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/correction.dat ]]; then
				echo -e "\texpected distribution calculation"
				/home/congm1/savanna/savannacong33/Code/SAD/bin/readsalmonbias correction ${OutDirectory}/salmon/${t}_${ID}_${n}/aux_info/ ${OutDirectory}/simu_${t}_${ID}_${n}_reference.fa ${OutDirectory}/salmon/${t}_${ID}_${n}/quant.sf ${OutDirectory}/salmon/${t}_${ID}_${n}/correction.dat
			fi

			# observed distribution
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/startpos.dat ]]; then
				echo -e "\tobserved distribution calculation"
				/home/congm1/savanna/savannacong33/Code/SAD/bin/transcovdist2 0 ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/salmon/${t}_${ID}_${n}/quant.sf ${OutDirectory}/salmon/${t}_${ID}_${n}/aux_info/eq_classes.txt ${OutDirectory}/salmon/${t}_${ID}_${n}/mapping.bam ${OutDirectory}/salmon/${t}_${ID}_${n}/startpos.dat
			fi

			# expression difference: comparing simulated reads info and salmon quant
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/expression_difference.txt ]]; then
				echo -e "\tcomparing simulated ground truth with salmon quant"
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/GroupTruthExp.py ${OutDirectory}/simu_${t}_${ID}_${n}_target.fa ${OutDirectory}/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${OutDirectory}/salmon/${t}_${ID}_${n}/quant.sf ${gtffile} ${OutDirectory}/salmon/${t}_${ID}_${n}/
			fi

			# run SAD
			if [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall ]]; then
				 echo -e "\tRunning SAD"
				mkdir -p ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/
				/home/congm1/savanna/savannacong33/Code/SAD/bin/SAD ${OutDirectory}/simu_${t}_${ID}_${n}_annotation.gtf ${OutDirectory}/salmon/${t}_${ID}_${n}/quant.sf ${OutDirectory}/salmon/${t}_${ID}_${n}/correction.dat ${OutDirectory}/salmon/${t}_${ID}_${n}/startpos.dat ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test
			fi

			# post-process SAD
			if [[ -e ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall ]] && [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted ]]; then
				echo -e "\tpost-process SAD p values"
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/AdjustPValue.py 1 ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted
			fi

			# evaluate SAD new isoform prediction
			if [[ -e ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted ]] && [[ ! -e ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_overall_exist ]]; then
				echo -e "\tevaluating SAD new isoform prediction"
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted_uniqgene
				python3 /home/congm1/savanna/savannacong33/Code/SAD/src/Simulation.py evatrans ${OutDirectory}/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${OutDirectory}/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/test_pvalue_overall_sorted_uniqgene ${gtffile} ${OutDirectory}/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_overall_exist
			fi

		done
	done
done
