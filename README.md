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
After obtaining the prerequisite libraries, SAD can be compiled by
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


## Specification of individual executable
### Retrieve the expected coverage distribution
### Retrieve the observed coverage distribution
### Detect and categorize anomalies
