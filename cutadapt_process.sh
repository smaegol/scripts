#!/bin/bash

## Preprocessing of fastq reads from NextSeq using cutadapt. Trims adapters, flter by quality and length.

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

### other variables
##cutadapt options
#adapters sequences
#TruSeqHT (double indexing)
adapter1='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT'  #for R1 reads 
adapter2='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG'  #for R2 reads 

cutadapt_error_rate='0.1' #error rate in overlapping fragment
cutadapt_count='2' #number of adapters
cutadapt_overlap='3' #length of min adapter overlap

min_quality='20' #quality threshold 
min_length='30' #minimum length after trimming

threads='20' #number of parallel threads to use 



#Get variables from command line
while getopts ":i:t:h:" optname
  do
    case "$optname" in
	"h")
		echo "Cutadapt processing"
		echo ""
		echo "Options :"
		echo " -i : input directory"
		echo " -t : threads"
		echo " -h : THIS HELP"
		echo ""
		exit 1
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


#### load modules (on the server I'm using environmental modules. Will load the most recent version of cutadapt (tested on 1.18)
module load cutadapt

for R1_FILE in `find $INPUT_DIR -name "*R1*.fastq*"`
do
	
	#R1_FILENAME=`expr match "$R1_FILE" '.*\/\(.*\)'`
	R1_PREFIX=`expr match "$R1_FILE" '\(.*\)R1.*'`
	R1_SUFFIX=`expr match "$R1_FILE" '.*R1\(.*\)'`
	R2_FILE=$R1_PREFIX"R2"$R1_SUFFIX
	R1_CUTADAPT_OUTPUT=$R1_PREFIX"R1_cutadapt_q"$min_quality"_l"$min_length".fastq.gz"
	R2_CUTADAPT_OUTPUT=$R1_PREFIX"R2_cutadapt_q"$min_quality"_l"$min_length".fastq.gz"
	cutadapt_log=$R1_PREFIX"cutadapt.log"
	echo "processing files: $R1_FILE and $R2_FILE"

	cutadapt -j $threads -a $adapter1 -A $adapter2 -o $R1_CUTADAPT_OUTPUT -p $R2_CUTADAPT_OUTPUT --nextseq-trim=$min_quality -m $min_length -n $cutadapt_count -e $cutadapt_error_rate -O $cutadapt_overlap $R1_FILE $R2_FILE &> $cutadapt_log

done
