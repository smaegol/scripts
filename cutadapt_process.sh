#!/bin/bash

## Preprocessing of fastq reads from NextSeq using cutadapt. Trims adapters, flter by quality and length. Process all fastq[.gz] files in given folder (and subfolders)

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
OUTPUT_DIR="." # default output dir - current dir
LIBRARY_TYPE="TruSeqHT" # default library type
find_maxdepth=1
filter_quality=1
#Arrays containing adapter sequences
declare -A adapter2=( 
 [TruSeqHT]="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"  
 [Nextera]="CTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT"  
 [dUTP]="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
 [TruSeqUD]="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
)
declare -A adapter1=( 
 [TruSeqHT]="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"  
 [Nextera]="CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"  
 [dUTP]="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
 [TruSeqUD]="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
)

### other variables
cutadapt_error_rate='0.1' #error rate in overlapping fragment
cutadapt_count='3' #number of adapters
cutadapt_overlap='3' #length of min adapter overlap

min_quality='20' #quality threshold 
min_length='30' #minimum length after trimming

threads='20' #number of parallel threads to use 

NextSeq=0
paired=1

function help {
	echo "Cutadapt processing of paired fastq files in given folder"
	echo ""
	echo "Options :"
	echo " -i : input directory [required]. Can be set as . to process current dir"
	echo " -o : output dir [defaults to input dir]"
	echo " -t : threads [default=$threads]"
	echo " -l : library type (one of TruSeqHT, Nextera, dUTP, TruSeqUD) [default=$LIBRARY_TYPE]"
	echo " -p : are data paired? (default=1)"
	echo " -q : min quality for filtering [default=$min_quality]"
	echo " -m : min length (-m option in cutadapt) [default=$min_length]"
	echo " -n : adapter count (-n option in cutadapt) [default=$cutadapt_count]"
	echo " -N : are data coming from NextSeq (default=$NextSeq, force use of --nextseq_trim instead of -q)"
	echo " -d : -maxdepth option for find - how deep go into folder structure [default = $find_maxdepth]" 
	echo " -f : filer by quality and length [default = $filter_quality, set to 0 to trim adapters only]"
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":hi:o:l:t:q:n:m:d:p:f:N:" optname
  do
    case "$optname" in
	"h")
		help
		exit 1
		;;
	"i")
		INPUT_DIR=$OPTARG
        	;;
	"o")
		OUTPUT_DIR=$OPTARG
        	;;
	"l")
		LIBRARY_TYPE=$OPTARG
        	;;
	"t")
        	threads=$OPTARG
	        ;;
	"q")
        	min_quality=$OPTARG
	        ;;
	"n")
        	cutadapt_count=$OPTARG
	        ;;
	"m")
        	min_length=$OPTARG
	        ;;
	"d")
        	find_maxdepth=$OPTARG
	        ;;
	"f")
        	filter_quality=$OPTARG
	        ;;
	"p")
        	paired=$OPTARG
	        ;;
	"N")
        	NextSeq=$OPTARG
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
		exit 1
		;;
   esac
done
if [ $OPTIND -eq 1 ]; then 
	echo "Please provide options"
	help
	exit 1 
fi

if [ -n "${adapter1[$LIBRARY_TYPE]}" ]; then
	adapter1=${adapter1[$LIBRARY_TYPE]}
	adapter2=${adapter2[$LIBRARY_TYPE]}
else
	echo "Wrong library type provided. Exiting..."
	exit 1
fi

# if using default OUTPUT DIR - set it to INPUT DIR
if [ "$OUTPUT_DIR" == "." ]; then
	OUTPUT_DIR=$INPUT_DIR
fi

if [ ! -d "$OUTPUT_DIR" ]; then
	echo "Output dir $OUTPUT_DIR does not exists. Creating now..."
	mkdir $OUTPUT_DIR
fi

if [ $NextSeq -eq 1 ]; then
	echo "Setting quality filtering for NextSeq data (--nextseq_trim==$min_quality)"
	quality_option="--nextseq-trim"
else
	echo "Setting quality filtering for data other than NextSeq (--quality-cutoff==$min_quality)"
	quality_option="--quality-cutoff"
fi



#### load modules (on the server I'm using environmental modules. Will load the most recent version of cutadapt (tested on 1.18)
module load cutadapt
module load fastqc

for R1_FILE in `find $INPUT_DIR -maxdepth $find_maxdepth -name "*R1*.fastq*"`
do
	R1_FILENAME=`expr match "$R1_FILE" '.*\/\(.*\)'`
	R1_PREFIX=`expr match "$R1_FILENAME" '\(.*\)R1.*'`
	R1_SUFFIX=`expr match "$R1_FILENAME" '.*R1\(.*\)'`
	if [ $paired -eq 1 ]; then
		R2_FILE=$INPUT_DIR"/"$R1_PREFIX"R2"$R1_SUFFIX
	fi
	if [ $filter_quality -eq 1 ]; then
		R1_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R1_cutadapt_q"$min_quality"_l"$min_length".fastq.gz"
		if [ $paired -eq 1 ]; then
			R2_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R2_cutadapt_q"$min_quality"_l"$min_length".fastq.gz"
		fi
	else 
		R1_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R1_cutadapt.fastq.gz"
		if [ $paired -eq 1 ]; then
			R2_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R2_cutadapt.fastq.gz"
		fi
	fi
	cutadapt_log=$OUTPUT_DIR"/"$R1_PREFIX"cutadapt.log"
	if [ $paired -eq 1 ]; then
		echo "processing files: $R1_FILE and $R2_FILE"
	else 
		echo "processing file: $R1_FILE"
	fi

	if [ $filter_quality -eq 1 ]; then
		if [ $paired -eq 1 ]; then
			if [ ! -e "$R2_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
				cutadapt -j $threads -a $adapter1 -A $adapter2 -o $R1_CUTADAPT_OUTPUT -p $R2_CUTADAPT_OUTPUT $quality_option=$min_quality -m $min_length -n $cutadapt_count -e $cutadapt_error_rate -O $cutadapt_overlap $R1_FILE $R2_FILE &> $cutadapt_log
			fi
		else
			if [ ! -e "$R1_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
				cutadapt -j $threads -a $adapter1 -o $R1_CUTADAPT_OUTPUT $quality_option=$min_quality -m $min_length -n $cutadapt_count -e $cutadapt_error_rate -O $cutadapt_overlap $R1_FILE &> $cutadapt_log
			fi
		fi
	else
		if [ $paired -eq 1 ]; then
			if [ ! -e "$R2_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
				cutadapt -j $threads -a $adapter1 -A $adapter2 -o $R1_CUTADAPT_OUTPUT -p $R2_CUTADAPT_OUTPUT -n $cutadapt_count -e $cutadapt_error_rate -O $cutadapt_overlap $R1_FILE $R2_FILE &> $cutadapt_log
			fi
		else
			if [ ! -e "$R1_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
				cutadapt -j $threads -a $adapter1 -o $R1_CUTADAPT_OUTPUT -n $cutadapt_count -e $cutadapt_error_rate -O $cutadapt_overlap $R1_FILE  &> $cutadapt_log
			fi
		fi
	fi
done
#Perform quality check after filtering
fastqc -t $threads --noextract $OUTPUT_DIR/*cutadapt*.fastq.gz
