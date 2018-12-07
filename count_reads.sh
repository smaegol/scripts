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
GENCODE_REF="/home/smaegol/data/indexes/Mus/annotation/gencode.vM19.annotation.gtf"

### other variables
format="bam"
order="pos"
feature_type="exon"
id_attribute="gene_id"
additional_attribute="gene_name"
mode="union"
nonunique="none"
stranded="yes"
threads='20' #number of parallel threads to use 


function help {
	echo "Count reads using HTSeq-count"
	echo ""
	echo "Options :"
	echo " -i : input directory [required]. Can be set as . to process current dir"
	echo " -o : output dir [defaults to input dir]"
	echo " -t : threads [default=$threads]"
	echo " -s : stranded [yes/no/reverse]"
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
while getopts ":hi:o:l:t:q:n:m:d:f:N:" optname
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

# if using default OUTPUT DIR - set it to INPUT DIR
if [ "$OUTPUT_DIR" == "." ]; then
	OUTPUT_DIR=$INPUT_DIR
fi

if [ ! -d "$OUTPUT_DIR" ]; then
	echo "Output dir $OUTPUT_DIR does not exists. Creating now..."
	mkdir $OUTPUT_DIR
fi

LOG_FILE=$OUTPUT_DIR"/count_reads.log"


#### load modules (on the server I'm using environmental modules. Will load the most recent version of cutadapt (tested on 1.18)
module load htseq

for BAM_FILE in `find $INPUT_DIR -name "*.bam"`
do
	BAM_FILENAME=`expr match "$BAM_FILE" '.*\/\(.*\)'`
	BAM_PREFIX=`expr match "$BAM_FILENAME" '\(.*\).bam'`
	HTSEQ_OUT=$OUTPUT_DIR"/"$BAM_PREFIX".htseq.txt"
	echo "$BAM_FILE"
	htseq-count -m $mode -f $format -r $order -s $stranded -t $feature_type -i $id_attribute --additional-attr $additional_attribute --nonunique $nonunique $BAM_FILE $GENCODE_REF > $HTSEQ_OUT
	
done
#Perform quality check after filtering

