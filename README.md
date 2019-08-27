# Overview
Salmon Anomaly Detection (SAD) detects the potential misquantifications for the RNA-seq transcript expression estimation made by Salmon. SAD detects large deviation of the observed coverage distribution from the expected coverage distribution for each transcript, and use the deviation as an indicator of misquantifications.

# Installation
## Prerequisite
SAD depends on the following libraries:
+ [Boost](https://www.boost.org/)
+ [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ [GSL](https://www.gnu.org/software/gsl/)
+ [Jellyfish](https://github.com/gmarcais/Jellyfish)
+ [HTSlib](http://www.htslib.org/)
+ [Spline](https://kluge.in-chemnitz.de/opensource/spline/)
+ [GUROBI](http://www.gurobi.com/)

For linux machine, the following script can be used to download and install the above dependencies EXCEPT GUROBI in current directory. Since GUROBI requires either commercial or academic license, we do not provide script for installation.
```
./install.sh
```

## Compiling SAD
After obtaining the prerequisite libraries, update pkg-config path to include the directory where the pkg-config files of dependencies are installed. If you install the dependencies using `install.sh` script, the pkg-config path can be updated by running the following command in shell:
```
export PKG_CONFIG_PATH=<path to the installation directory>/external/pkgconfig:${PKG_CONFIG_PATH}
```
SAD can be compiled by
```
cd sad/
make
```

# Usage
Three steps are needed for SAD: retrieve the observed coverage distribution, retrieve the expected coverage distribution, detect and categorize anomalies.

We provide a python script to go through all the steps. Alternatively, three executables are provided for the three steps, and users can run each individual executables.

## Quick start
The following python script will sequentially run the required steps, and generate the predicted unadjustable anomalies and the adjusted quantifications.
```
python3 SADpipe.py -t <transcriptome.fa> -a <annotation.gtf> -s <salmon_folder> -o <output_prefix>
``` 

The input Salmon quantification result should be generated with the following options:
```
salmon quant -i <salmon index> -1 <read_1> (-2 <reads_2>) --gcBias --seqBias --posBias --dumpEqWeights -o <output folder> --writeMappings=<output folder>/mapping.sam
```
With adding the above options, Salmon output directory should have the following structure, with which users can check their output:
- quant.sf
- mapping.sam
- aux_info/
	- fld.gz
	- exp_gc.gz
	- obs_gc.gz
	- exp5_pos.gz
	- obs5_pos.gz
	- exp3_pos.gz
	- obs3_pos.gz
	- exp5_seq.gz
	- obs5_seq.gz
	- exp3_seq.gz
	- obs3_seq.gz

## SAD output
Two files are the main outputs of SAD, corresponding to the p-value of unadjustable anomalies, and the adjusted quantification of the adjustable anomalies.
+ output_prefix_unadjustable_pvalue.tsv: p-value of unadjustable anomalies. Note that the p-values of all evaluated transcripts within the output, and users can define their own p-value cutoff.
	1. Name: the ID of the transcript.
	2. Coverage: defined as number reads / transcript length.
	3. AnomalyScoreNeg: transcript-level under-expression anomaly score.
	4. RegionStartNeg: the starting position of under-expression region in transcript coordinate.
	5. RegionEndNeg: the ending position of under-expression region in transcript coordinate.
	6. PValue_Neg: p-value of the transcript-level under-expression anomaly score
	7. AdjPValue_Neg: FDR adjusted p-value of the transcript-level under-expression anomaly score
	8. AnomalyScorePos: transcript-level over-expression anomaly score.
	9. RegionStartPos: the starting position of over-expression region in transcript coordinate.
	10. RegionEndPos: the ending position of over-expression region in transcript coordinate.
	11. PValue_Pos: p-value of the transcript-level over-expression anomaly score.
	12. AdjPValue_Pos: FDR adjusted p-value of the transcript-level over-expression anomaly score.
	13. MinAdjPValue: the minimum between the FDR adjusted under-expression anomaly p-value and over-expression anomaly p-value.
	14. Choice: indicator of which of the under-expression p-value and over-expression p-value is smaller. 0 refers to under-expression, 1 refers to over-expression.

+ output_prefix_adjusted_quantification.tsv: the adjusted quantification of the adjustable anomalies
	1. Name: the ID the transcript.
	2. Length: the length of the transcript.
	3. NumReads: the weighted number of reads assigned to the transcript by SAD re-assignment procedure.

## Specification of individual executable
### Retrieve the expected coverage distribution
After compiling SAD, bin/readsalmonbias is the executable of retrieving the expected distribution for each transcript using Salmon bias correction model. The `--gcBias`, `--seqBias`, `--posBias` options in Salmon infer the corresponding bias models. bin/readrsembias is the executable of retrieving the expected distribution under RSEM bias correction model. The `--estimate-rspd` option in RSEM allows it to estimate a positional bias. For both executables, the order of the input arguments matters. To run the executables, using
```
bin/readsalmonbias correction <salmon aux folder path> <transcriptome.fa> <salmon quant.sf> <output file> (number_threads)
```
or for RSEM,  locate the .model file under .stat folder and run
```
bin/readrsembias correction <rsem .model file> <transcriptome.fa> <output file>
```

The output is a binary file of the below content. If you want to use SAD for other quantifiers, outputing an expected coverage distribution file by a custom script is needed.

	- the number of transcripts (int32_t)
	- for each transcript:
		- the length of the transcript ID (int32_t)
		- transcript sequence length (int32_t)
		- transcript ID (char * length of ID)
		- expected coverage distribution (double * length of sequence)

### Retrieve the observed coverage distribution
bin/transcovdist is the executable of retrieving from Salmon the observed distribution for each transcript using the mapping-to-transcriptome BAM file. With --writeMappings option, Salmon will output a SAM file of the read mapping, and can be converted to BAM using samtools. bin/rsemobs is the exectuable of retrieving from RSEM the observed distribution. For both executables, the order of the input arguments matters. To run the executables, using
```
bin/transcovdist <annotation gtf file> <salmon quant.sf> <salmon eq_classes.txt> <salmon mapping.bam> <output file> 
```
or
```
bin/rsemobs <annotation gtf file> <RSEM mapping bam file> <output file>
```

The output is a binary file of below content. For other quantifiers, a custom script is required to output the observed coverage in the following structure in order to apply SAD.

	- the number of transcripts (int32_t)
	- for each transcript:
		- the length of the transcript ID (int32_t)
		- the number of effective positions (int32_t). Effetive positions are the positions with non-zero fragment starts besides the last position of transcript.
		- transcript ID (char * length of ID)
		- effective positions (int32_t * number of effective positions)
		- weighted count of fragment starts at each effective position (double * number of effective positions)

### Detect and categorize anomalies
bin/SAD is the main method for anomaly detection. It can be executed with the following command:
```
bin/SAD <annotation.gtf> <salmon quant.sf> <expected distribution> <observed distribution> <output prefix> (number_threads)
```
Note that the order of the input arguments matters.Please refer to the SAD output for the output specification for this executable.
