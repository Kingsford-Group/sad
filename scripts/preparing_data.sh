#!/bin/bash

folder=$0
folder=${folder%/*}

# download reference and build index
if [[ ! -e gencode.v26.annotation.gtf ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > gencode.v26.annotation.gtf
fi
if [[ ! -e GRCh38.p10.genome.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > GRCh38.p10.genome.fa
fi
if [[ ! -e gencode.v26.transcripts.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > gencode.v26.transcripts.fa
fi
if [[ ! -e gencode.v26.pc.transcripts.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.pc_transcripts.fa.gz | gunzip > gencode.v26.pc.transcripts.fa
fi
if [[ ! -e StarIndex/Genome ]]; then
	mkdir -p StarIndex
	STAR --runThreadN 8 --runMode genomeGenerate --genomeDir StarIndex/ --genomeFastaFiles GRCh38.p10.genome.fa
fi
if [[ ! -e gencode.v26.full/hash.bin ]]; then
	salmon index --gencode -t gencode.v26.transcripts.fa -i gencode.v26.full
fi
if [[ ! -e gencode.v26.pc/hash.bin ]]; then
	salmon index --gencode -t gencode.v26.pc.transcripts.fa -i gencode.v26.pc
fi
if [[ ! -e gencode.v26.Gene_Trans_Map.txt ]]; then
	awk '{if($3=="transcript") print substr($10,2,length($10)-3)"\t"substr($12,2,length($12)-3)}' gencode.v26.annotation.gtf > gencode.v26.Gene_Trans_Map.txt
fi

# extract protein-coding only transcript
python ${folder}/ExtractPCAnnotation.py gencode.v26.pc.transcripts.fa gencode.v26.annotation.gtf gencode.v26.annotation.pc.gtf

if false; then
# download data
echo "DOWNLOADING GEUVADIS SAMPLES..."
${folder}/../data/GEUVADIS/getdata.sh
echo "DOWNLOADING HUMAN BODY MAP SAMPLES..."
${folder}/../data/HumanBodyMap/getdata.sh
fi

# simulating deletion, fusion events of transcriptome; simulating reads with polyester
Type=("PC" "Full")
GTFfiles=("gencode.v26.annotation.pc.gtf" "gencode.v26.annotation.gtf")
TransFastafiles=("gencode.v26.pc.transcripts.fa" "gencode.v26.transcripts.fa")
NumEvents=(200 500 1000 1500)
RefQuantPrefix=("${folder}/../data/BaseExpression/ERR188297" "${folder}/../data/BaseExpression/SRR3192396" "${folder}/../data/BaseExpression/SRR3192412")
RefQuantID=(GEU GM12878 K562)
for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastafiles[${i}]}
	for n in ${NumEvents[@]}; do
		for ((j=0; j<${#RefQuantPrefix[@]}; j++)); do
			ID=${RefQuantID[${j}]}
			if ((i==0)); then
				refquantfile=${RefQuantPrefix[${j}]}"_salmon_pc_quant.sf"
			else
				refquantfile=${RefQuantPrefix[${j}]}"_salmon_full_quant.sf"
			fi

			echo simu_${t}_${ID}_${n}

			# simulating deletion and fusion
			if [[ ! -e simu_${t}_${ID}_${n}_theoexp.txt ]]; then
				echo -e "\tSimulating deletion and fusion."
				mkdir -p simu_${t}_${ID}_${n}
				python3 ${folder}/Simulation.py ${transfasta} ${gtffile} ${refquantfile} simu_${t}_${ID}_${n} ${n} ${n}
			fi

			# categorize the simulated events into with novel junctions and without novel junctions
			if [[ ! -e groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ]]; then
				echo -e "\tCategorizing simulated events."
				mkdir -p groundtruth_${t}_${ID}_${n}
				${folder}/../bin/categorizesimulation ${gtffile} simu_${t}_${ID}_${n}_removelist.txt simu_${t}_${ID}_${n}_target.fa groundtruth_${t}_${ID}_${n}/
			fi

			# simulating reads with polyester
			if [[ ! -e polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ]] && [[ ! -e  polyester_${t}_${ID}_${n}/sample_01_1.fasta ]]; then
				echo -e "\tSimulating reads with polyester."
				mkdir -p polyester_${t}_${ID}_${n}
				#/opt/local/stow/R-3.3.2/bin/Rscript /mnt/disk33/user/congm1/Code/SAD/src/SimulationReads.R simu_${t}_${ID}_${n}_target.fa simu_${t}_${ID}_${n}_theoexp.txt polyester_${t}_${ID}_${n}/ 1
				Rscript /mnt/disk33/user/congm1/Code/SAD/src/SimulationReads.R simu_${t}_${ID}_${n}_target.fa simu_${t}_${ID}_${n}_theoexp.txt polyester_${t}_${ID}_${n}/ 1
			fi

			# gzipping simulated reads
			if [[ $(find polyester_${t}_${ID}_${n}/ -name *fasta) != "" ]]; then
				echo -e "\tGzipping simulated reads."
				gzip polyester_${t}_${ID}_${n}/*fasta
			fi

		done
	done
done
