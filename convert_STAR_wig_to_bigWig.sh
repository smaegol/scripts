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
CHROM_SIZES="/home/smaegol/data/indexes/Mus/chromSizes"
OUT_DIR=$INPUT_DIR #default output dir

### other variables

find_maxdepth=1

additional_suffix=""
name_search="*Signal.UniqueMultiple.str*out.wig" #how to find wig files
include_multimappers=1
additional_suffix=""


function help {
	echo "STAR mapping script"
	echo ""
	echo "Options :"
	echo " -i : Input dir [defaults to current dir]"
	echo " -o : Output dir [defaults to input dir]"
	echo " -c : chromSize, other than default [$CHROM_SIZES]"
	echo " -m : include multimappers in bw [default=$include_multimappers]"
	echo " -a : additional suffix to add to output files (to differentiate different runs) [default = empty]"
	echo " -d : -maxdepth option for find - how deep go into folder structure [default = $find_maxdepth]" 
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":ho:i:c:a:m:d:" optname
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
	"c")
        	CHROM_SIZES=$OPTARG
	        ;;
	"a")
        	additional_suffix=$OPTARG
	        ;;
	"m")
        	include_multimappers=$OPTARG
	        ;;
	"d")
        	find_maxdepth=$OPTARG
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

if [ "$OUT_DIR" == "." ]; then
	OUT_DIR=$INPUT_DIR
fi

if [ ! -d "$OUT_DIR" ]; then
	mkdir $OUT_DIR
fi


#if dont want multimappers output
if [ $include_multimappers -eq 0 ]; then
	name_search="*Signal.Unique.str*out.wig"
fi


#### load modules (on the server I'm using environmental modules. Will load the most recent version of software
module load samtools
module load ucsc_tools

#### main part
for wig_file in `find $INPUT_DIR -maxdepth $find_maxdepth -name "$name_search"`
do
	WIG_FILENAME=`expr match "$wig_file" '.*\/\(.*\)'`
	WIG_PREFIX=`expr match "$WIG_FILENAME" '\(.*\)_Signal.*'`
	WIG_SUFFIX=`expr match "$WIG_FILENAME" '.*Signal.\(.*\)'`
	WIG_STRAND=`expr match "$WIG_FILENAME" '.*str\(.\)\..*'`
	if [ "$WIG_STRAND" == '1' ] ; then
		strand_to_write="rev" #default for dUTP
	else
		strand_to_write="fwd"
	fi
	bw_output=$OUT_DIR"/"$WIG_PREFIX""$strand_to_write$additional_suffix".bw"
	if [ ! -e "$bw_output" ]; 
	then
		echo "converting $wig_file to $bw_output..."
		wigToBigWig $wig_file $CHROM_SIZES $bw_output
	fi
done
echo "Finished processing"
