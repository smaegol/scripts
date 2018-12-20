#!/bin/bash

## Quality control of mapped data

#######################################################################################
###                                                                                 ###
###     Copyright (C) 2018  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

## variables ##
### Global variables - do not modify
INPUT_DIR="." # default input dir - current dir
GENOME_DIR="/home/smaegol/data/indexes/Mus/star_indices_overhang100/"
GENE_BODY_REF="/home/smaegol/data/indexes/Mus/annotation/RSEQC/mm10_Gencode_VM18.bed"
GENCODE_REF="/home/smaegol/data/indexes/Mus/annotation/gencode.vM19.annotation.gtf"
HOUSEKEEPING_BODY_REF="/home/smaegol/data/indexes/Mus/annotation/RSEQC/mm10.HouseKeepingGenes.bed"
QORTS_DIR="QORTS_QC"
### other variables
##cutadapt options
#adapters sequences
#TruSeqHT (double indexing)

threads='20' #number of parallel threads to use 

function help {
	echo "Quality control functions (RSeQC, DeepTools. Samtools, QoRTs)"
	echo "Creates quality control data, plots in the INPUT_DIR)"
	echo ""
	echo "Options :"
	echo " -i : input directory"
	echo " -t : threads"
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":o:i:t:h:" optname
  do
    case "$optname" in
	"h")
		help
		exit 1
		;;
	"o")
		OUT_DIR=$OPTARG
		;;
	"i")
		INPUT_DIR=$OPTARG
        	;;
	"t")
        	threads=$OPTARG
	        ;;
	"?")
        	echo "Unknown option $OPTARG"
        	exit 1 
        	;;
	":")
		echo "No argument value for option $OPTARG"
		exit 1
	        ;;
      *)
      # Should not occur
		echo "Unknown error while processing options"
		;;
   esac
done
if [ $OPTIND -eq 1 ]; then 
	echo "Please provide options"
	help
	exit 1 
fi

QORTS_OUTDIR=$INPUT_DIR"/"$QORTS_DIR
if [ ! -d "$QORTS_OUTDIR" ]; then
	mkdir "$QORTS_OUTDIR"
fi

#### load modules (on the server I'm using environmental modules. Will load the most recent version of software
module load samtools
module load deeptools
module load rseqc
module load R
module load qorts
#### main part

BAM_FILES=$(find $INPUT_DIR -name "*.bam"  | paste -d, -s)
BAM_FILES_COV=$(find $INPUT_DIR -name "*.bam"  | paste -d' ' -s)
PLOT_COV_OUT_RAW="deeptools_plot_coverage_raw_counts.txt"
PLOT_COV_OUT_PLOT="deeptools_plot_coverage.png"
GENE_BODY_COV_PREFIX="gene_body_cov"



for BAM_FILE in `find $INPUT_DIR -name "*.bam"`
do
	BAM_FILENAME=`expr match "$BAM_FILE" '.*\/\(.*\)'`
	BAM_PREFIX=`expr match "$BAM_FILENAME" '\(.*\).bam'`
	QORTS_OUT=$QORTS_OUTDIR"/"$BAM_PREFIX
	if [ ! -e "$BAM_FILE.bai" ] ; then
		echo "Bam file $BAM_FILE is not indexed. Generating index now..."
		samtools index $BAM_FILE
	fi
	#generate idxstats
	if [ ! -e $BAM_FILE".idxstats" ]; then
		echo "generating idxstats of $BAM_FILE"
		samtools idxstats -@ $threads $BAM_FILE > $BAM_FILE".idxstats"
	fi
	#RSEQC stats:
	if [ ! -e $BAM_FILE".infer_experiment" ]; then
		echo "inferring experiment of $BAM_FILE"
		infer_experiment.py -i $BAM_FILE -r $GENE_BODY_REF > $BAM_FILE".infer_experiment"
	fi
	if [ ! -e $BAM_FILE".read_distribution.txt" ]; then
		echo "calculating read distribution of $BAM_FILE"
		read_distribution.py -i $BAM_FILE -r $GENE_BODY_REF > $BAM_FILE".read_distribution.txt"
	fi
	if [ ! -e $BAM_FILE".inner_distance.txt" ]; then
		echo "summarizing inner distance of $BAM_FILE"
		inner_distance.py -i $BAM_FILE -r $GENE_BODY_REF  -o $BAM_FILE
	fi
	if [ ! -e $BAM_FILE".bam_stat.txt" ]; then
		echo "calculating bam_stat of $BAM_FILE"
		bam_stat.py -i $BAM_FILE > $BAM_FILE".bam_stat.txt"
	fi
	if [ ! -e $BAM_FILE".clipping_profile.r" ]; then
		echo "calculating clipping profile of $BAM_FILE"
		clipping_profile.py -i $BAM_FILE -o $BAM_FILE -s PE
	fi
	if [ ! -e $BAM_FILE".junctionSaturation_plot.r" ]; then
		echo "calculating junction saturation of $BAM_FILE"
		junction_saturation.py -i $BAM_FILE -r $GENE_BODY_REF -o $BAM_FILE
	fi
	if [ ! -d "$QORTS_OUT" ]; then
		mkdir "$QORTS_OUT"
		QoRTs.jar QC --generatePlots $BAM_FILE $GENCODE_REF $QORTS_OUT
	fi
done

if [ ! -e "$PLOT_COV_OUT_PLOT" ]; then 
	echo "calculating coverage of all files"
	plotCoverage --bamfiles $BAM_FILES_COV --outRawCounts $PLOT_COV_OUT_RAW -o $PLOT_COV_OUT_PLOT -p $threads
fi

if [ ! -e "$GENE_BODY_COV_PREFIX"".geneBodyCoverage.txt" ]; then
	echo "calculating gene body coverage of all files"
	geneBody_coverage.py -i $BAM_FILES -r $HOUSEKEEPING_BODY_REF -f png -o $GENE_BODY_COV_PREFIX
fi

if [ ! -e "tin_report.txt" ]; then
	echo "calculating tin of all files"
	tin.py -i $BAM_FILES -r $HOUSEKEEPING_BODY_REF > "tin_report.txt"
fi
