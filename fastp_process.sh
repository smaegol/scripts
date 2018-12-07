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
find_maxdepth=1
filter_quality=1
#Arrays containing adapter sequences

### other variables

min_quality='20' #quality threshold 
min_length='30' #minimum length after trimming

threads='20' #number of parallel threads to use 


function help {
	echo "Cutadapt processing of paired fastq files in given folder"
	echo ""
	echo "Options :"
	echo " -i : input directory [required]. Can be set as . to process current dir"
	echo " -o : output dir [defaults to input dir]"
	echo " -t : threads [default=$threads]"
	echo " -q : min quality for filtering [default=$min_quality]"
	echo " -m : min length (-m option in cutadapt) [default=$min_length]"
	echo " -d : -maxdepth option for find - how deep go into folder structure [default = $find_maxdepth]" 
	echo " -f : filer by quality and length [default = $filter_quality, set to 0 to trim adapters only]"
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":h:i:o:l:t:q:n:m:d:f:N:" optname
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
	"t")
        	threads=$OPTARG
	        ;;
	"q")
        	min_quality=$OPTARG
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




#### load modules (on the server I'm using environmental modules. Will load the most recent version of cutadapt (tested on 1.18)
module load cutadapt
module load fastqc
module load fastp

for R1_FILE in `find $INPUT_DIR -maxdepth $find_maxdepth -name "*R1*.fastq*"`
do
	R1_FILENAME=`expr match "$R1_FILE" '.*\/\(.*\)'`
	R1_PREFIX=`expr match "$R1_FILENAME" '\(.*\)R1.*'`
	R1_SUFFIX=`expr match "$R1_FILENAME" '.*R1\(.*\)'`
	R2_FILE=$INPUT_DIR"/"$R1_PREFIX"R2"$R1_SUFFIX
	if [ $filter_quality -eq 1 ]; then
		R1_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R1_fastp_q"$min_quality"_l"$min_length".fastq.gz"
		R2_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R2_fastp_q"$min_quality"_l"$min_length".fastq.gz"
	else 
		R1_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R1_fastp.fastq.gz"
		R2_CUTADAPT_OUTPUT=$OUTPUT_DIR"/"$R1_PREFIX"R2_fastp.fastq.gz"
	fi
	fastp_json=$OUTPUT_DIR"/"$R1_PREFIX"fastp.json"
	fastp_html=$OUTPUT_DIR"/"$R1_PREFIX"fastp.html"
	echo "processing files: $R1_FILE and $R2_FILE"

	if [ $filter_quality -eq 1 ]; then
		if [ ! -e "$R2_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
			fastp -w $threads -o $R1_CUTADAPT_OUTPUT -O $R2_CUTADAPT_OUTPUT -q $min_quality -l $min_length -i $R1_FILE -I $R2_FILE --cut_by_quality3 -p -c -j $fastp_json -h $fastp_html --detect_adapter_for_pe
		fi
	else
		if [ ! -e "$R2_CUTADAPT_OUTPUT" ]; then # Run cutadapt only if output is missing	
			fastp -w $threads -o $R1_CUTADAPT_OUTPUT -O $R2_CUTADAPT_OUTPUT  -i $R1_FILE -I $R2_FILE -p -j $fastp_json -h $fastp_html -Q -L --detect_adapter_for_pe
		fi
	fi
done
#Perform quality check after filtering
fastqc -t $threads --noextract $OUTPUT_DIR/*fastp*.fastq.gz
