#!/bin/bash

import sys
import subprocess
from pathlib import Path

# input:
# 	- transcript fasta file
# 	- gene annotation GTF file
# 	- salmon folder:
# 		- quant.sf
# 		- mapping.sam
# 		- aux/
# 	- output folder


def ParseArgument(arguments):
	Mode = 0
	TranscriptFasta = ""
	GTFfile = ""
	SalmonDir = ""
	RSEMPrefix = ""
	OutPrefix = ""
	i = 1
	while i < len(arguments):
		if arguments[i] == "-m":
			assert( i+1 < len(arguments) )
			if arguments[i+1] == "0":
				Mode = 0
			elif arguments[i+1] == "1":
				Mode = 1
		elif arguments[i] == "-t":
			assert( i+1 < len(arguments) )
			TranscriptFasta = arguments[i+1]
			print("TranscriptFasta = " + TranscriptFasta)
		elif arguments[i] == "-a":
			assert( i+1 < len(arguments) )
			GTFfile = arguments[i+1]
			print("GTFfile = " + GTFfile)
		elif arguments[i] == "-s":
			assert( i+1 < len(arguments) )
			SalmonDir = arguments[i+1]
			print("SalmonDir = " + SalmonDir)
		elif arguments[i] == "-o":
			assert( i+1 < len(arguments) )
			OutPrefix = arguments[i+1]
			print("OutPrefix = " + OutPrefix)
		else:
			print("Error: Invalid argument {}".format(arguments[i]))
			sys.exit()
		i += 2
	if TranscriptFasta == "":
		print("Please input transcriptome fasta file.")
		sys.exit()
	if GTFfile == "":
		print("Please input annotation GTF file.")
		sys.exit()
	if Mode == 0 and SalmonDir == "":
		print("Please input the directory of Salmon / RSEM quantification.")
		sys.exit()
	if Mode == 1 and RSEMPrefix == "":
		print("Please input the RSEM output prefix.")
		sys.exit()
	if OutPrefix == "":
		print("Please input the output prefix")
		sys.exit()
	return Mode, TranscriptFasta, GTFfile, SalmonDir, RSEMPrefix, OutPrefix


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python3 SADpipe.py -m <mode (0 or 1)> -t <transcriptome.fa> -a <annotation.gtf> -s <salmon_folder when mode 0 or rsem_prefix when mode 1> -o <output_prefix>")
		print("\tmode 0: for Salmon quantifier (default)")
		print("\tmode 1: for RSEM quantifier")
	else:
		codedir = "/".join(sys.argv[0].split("/")[:-1])
		if codedir == "":
			codedir = "."
		Mode, TranscriptFasta, GTFfile, SalmonDir, RSEMPrefix, OutPrefix = ParseArgument(sys.argv)

		if Mode == 0:
			# process bias correction
			assert( Path(SalmonDir + "/aux_info/fld.gz").exists() )
			assert( Path(SalmonDir + "/aux_info/exp_gc.gz").exists() and Path(SalmonDir + "/aux_info/obs_gc.gz").exists() )
			assert( Path(SalmonDir + "/aux_info/exp5_pos.gz").exists() and Path(SalmonDir + "/aux_info/obs5_pos.gz").exists() and Path(SalmonDir + "/aux_info/exp3_pos.gz").exists() and Path(SalmonDir + "/aux_info/obs3_pos.gz").exists() )
			assert( Path(SalmonDir + "/aux_info/exp5_seq.gz").exists() and Path(SalmonDir + "/aux_info/obs5_seq.gz").exists() and Path(SalmonDir + "/aux_info/exp3_seq.gz").exists() and Path(SalmonDir + "/aux_info/obs3_seq.gz").exists() )
			if not Path(SalmonDir + "/correction.dat").exists():
				print("PROCESSING EXPECTED COVERAGE DISTRIBUTION...")
				aux_dir = SalmonDir + "/aux_info"
				quant_file = SalmonDir + "/quant.sf"
				output_file = SalmonDir + "/correction.dat"
				p = subprocess.Popen("{}/../bin/readsalmonbias correction {} {} {} {}".format(codedir, aux_dir, TranscriptFasta, quant_file, output_file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()

			# converting SAM to BAM
			if Path(SalmonDir + "/mapping.sam").exists() and not Path(SalmonDir + "/mapping.bam").exists():
				print("CONVERTING SAM TO BAM...")
				p = subprocess.Popen("samtools view -Shb {} -o {}".format(SalmonDir+"/mapping.sam", SalmonDir+"/mapping.bam"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()

			# processing observed distribution
			assert( Path(SalmonDir + "/mapping.bam").exists() )
			if not Path(SalmonDir + "/startpos.dat").exists():
				print("PROCESSING OBSERVED COVERAGE DISTRIBUTION...")
				quant_file = SalmonDir + "/quant.sf"
				eq_file = SalmonDir + "/aux_info/eq_classes.txt"
				bam_file = SalmonDir + "/mapping.bam"
				output_file = SalmonDir + "/startpos.dat"
				p = subprocess.Popen("{}/../bin/transcovdist 0 {} {} {} {} {}".format(codedir, GTFfile, quant_file, eq_file, bam_file, output_file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()

			# running SAD
			if not Path(OutPrefix + "_unadjustable_pvalue.tsv").exists():
				print("DETECTING ANOMALIES USING SAD...")
				# creating output directory
				p = subprocess.Popen("mkdir -p {}".format("/".join(OutPrefix.split("/")[:-1])), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()
				# running SAD
				quant_file = SalmonDir + "/quant.sf"
				correction_file = SalmonDir + "/correction.dat"
				startpos_file = SalmonDir + "/startpos.dat"
				p = subprocess.Popen("{}/../bin/SAD 0 {} {} {} {} {}".format(codedir, GTFfile, quant_file, correction_file, startpos_file, OutPrefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()
		else:
			RSEMFolder = "/".join(RSEMPrefix.split("/")[:-1])
			# processing bias correction for RSEM
			if not Path(RSEMFolder + "/correction.dat").exists():
				print("PROCESSING EXPECTED COVERAGE DISTRIBUTION...")
				p = subprocess.Popen("{}/../bin/readrsembias correction {} {} {}".format(codedir, RSEMPrefix + ".model", TranscriptFasta, RSEMFolder + "/correction.dat"))
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()

			# processing observed distribution
			if not Path(RSEMFolder + "/startpos.dat").exists():
				print("PROCESSING OBSERVED COVERAGE DISTRIBUTION...")
				p = subprocess.Popen("{}/../bin/rsemobs {} {}".format(codedir, RSEMPrefix + ".transcript.bam", RSEMFolder + "/startpos.dat"))
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()

			# running SAD
			if not Path(OutPrefix + "_unadjustable_pvalue.tsv").exists():
				print("DETECTING ANOMALIES USING SAD...")
				# creating output directory
				p = subprocess.Popen("mkdir -p {}".format("/".join(OutPrefix.split("/")[:-1])), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()
				# running SAD
				quant_file = RSEMPrefix + ".isoforms.results"
				correction_file = RSEMFolder + "/correction.dat"
				startpos_file = RSEMFolder + "/startpos.dat"
				p = subprocess.Popen("{}/../bin/SAD 1 {} {} {} {} {}".format(codedir, GTFfile, quant_file, correction_file, startpos_file, OutPrefix), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate()
				if err != b'':
					print(err)
					sys.exit()