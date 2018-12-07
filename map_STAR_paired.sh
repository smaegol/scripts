#!/bin/bash

## Mapping or reads using STAR (single invovation for every reads pair)

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
OUT_DIR=$INPUT_DIR #default output dir
### other variables

threads='20' #number of parallel threads to use 

function help {
	echo "STAR mapping script"
	echo ""
	echo "Options :"
	echo " -1 : First input file"
	echo " -2 : Second input file"
	echo " -o : Output dir"
	echo " -t : threads"
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":1:2:o:i:t:h:" optname
  do
    case "$optname" in
	"h")
		help
		exit 1
		;;
	"1")
		R1_FILE=$OPTARG
		;;
	"2")
		R2_FILE=$OPTARG
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

if [ -z $R1_FILE ] || [ -z $R2_FILE ]; then
	echo "Input files missing, please provide -1 and -2 options"
	exit 1
fi

if [ ! -d "$OUT_DIR" ]; then
	mkdir $OUT_DIR
fi

#### load modules (on the server I'm using environmental modules. Will load the most recent version of software
module load star
module load samtools

#### main part

R1_PREFIX=`expr match "$R1_FILE" '\(.*\)R1.*'`
R1_SUFFIX=`expr match "$R1_FILE" '.*R1\(.*\)'`
STAR_PREFIX=$OUT_DIR"/"$R1_PREFIX
STAR_SAM_OUTPUT=$STAR_PREFIX"Aligned.out.sam"
STAR_BAM_OUTPUT=$STAR_PREFIX"Aligned.out.sorted.bam"
echo "Mapping files: $R1_FILE and $R2_FILE"
#map using STAR:
if [ ! -e "$STAR_SAM_OUTPUT" ]; 
then
	STAR --runThreadN $threads --genomeDir $GENOME_DIR --readFilesIn $R1_FILE $R2_FILE --readFilesCommand gunzip -c --outFileNamePrefix $STAR_PREFIX --quantMode GeneCounts --outFilterMultimapNmax 1000 
fi
#convert to sorted bam
if [ ! -e "$STAR_BAM_OUTPUT" ]; 
then
	echo "Converting to sorted bam..."
	samtools view -b $STAR_SAM_OUTPUT | samtools sort -@ $threads -o $STAR_BAM_OUTPUT
	samtools index $STAR_BAM_OUTPUT
fi
echo "Finished processing"
