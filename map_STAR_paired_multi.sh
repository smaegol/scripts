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
genome_load_option="LoadAndKeep" # force STAR to keep loaded genome in the memory till the end of all mapping tasks

### other variables
quant_mode="--quantMode GeneCounts" # handy for further multiqc
convert_to_sorted_bam=0
additional_star_options=""
gene_counts=1
threads='20' #number of parallel threads to use 
find_maxdepth=1
read_files_command="--readFilesCommand gunzip -c"
out_wig_type="wiggle"
out_wig_norm="RPM"
out_sam_mode="Full"
out_sam_type="BAM SortedByCoordinate"
out_sam_strand_field="intronMotif"
out_reads_unmapped="--outReadsUnmapped Fastx" #output unmapped reads? default yes
output_unmapped=1
limit_bam_sort_ram=5000000000 # required to be >0 if using genomeLoad other than NoShared and outputting Sorted Bam
output_wiggle=1 # default to 1 - outputs wiggle (take some time and disk space)
additional_suffix=""
paired=1

function help {
	echo "STAR mapping script"
	echo ""
	echo "Options :"
	echo " -i : Input dir [defaults to current dir]"
	echo " -o : Output dir [defaults to input dir]"
	echo " -g : Genome dir for STAR, other than default [$GENOME_DIR]"
	echo " -t : threads [default 20]"
	echo " -c : perform gene counts [default=1]"
	echo " -a : additional STAR options [enclose in \" \"]"
#	echo " -s : convert to sorted bam using samtools [default=1]"
	echo " -d : -maxdepth option for find - how deep go into folder structure [default = $find_maxdepth]" 
	echo " -p : additional suffix to add to output files (to differentiate different runs) [default = empty]"
	echo " -r : are reads paired? (default=1)"
	echo " -u : output unmapped reads [default = 1]"
	echo " -w : output wiggle [default = 1]"
	echo " -h : THIS HELP"
	echo ""
}

#Get variables from command line
while getopts ":o:i:t:c:a:g:hd:p:r:u:w:" optname
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
	"c")
        	gene_counts=$OPTARG
	        ;;
	"a")
        	additional_star_options=$OPTARG
	        ;;
	"g")
        	GENOME_DIR=$OPTARG
	        ;;
	"d")
        	find_maxdepth=$OPTARG
	        ;;
	"w")
        	output_wiggle=$OPTARG
	        ;;
	"p")
        	additional_suffix=$OPTARG
	        ;;
	"r")
        	paired=$OPTARG
	        ;;
	"u")
        	output_unmapped=$OPTARG
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

if [ $gene_counts -eq 0 ]; then
	quant_mode=""
fi

#if dont want unmapped output
if [ $output_unmapped -eq 0 ]; then
	out_reads_unmapped=""
fi

if [ $output_wiggle -eq 1 ]; then
	wiggle_output_command="--outWigType $out_wig_type --outWigNorm $out_wig_norm"
else
	wiggle_output_command=""
fi

#### load modules (on the server I'm using environmental modules. Will load the most recent version of software
module load star
module load samtools

# load the genome for mapping
STAR --genomeDir $GENOME_DIR --genomeLoad LoadAndExit
#### main part
for R1_FILE in `find $INPUT_DIR -maxdepth $find_maxdepth -name "*R1*.fastq*"`
do
	R1_FILENAME=`expr match "$R1_FILE" '.*\/\(.*\)'`	
	R1_PREFIX=`expr match "$R1_FILENAME" '\(.*\)R1.*'`
	R1_SUFFIX=`expr match "$R1_FILENAME" '.*R1\(.*\)'`
	if [ $paired -eq 1 ]; then
		R2_FILE=$INPUT_DIR"/"$R1_PREFIX"R2"$R1_SUFFIX
	else
		R2_FILE=''
	fi
	STAR_PREFIX=$OUT_DIR"/"$R1_PREFIX"_"$additional_suffix
	if [ "$out_sam_type" == "BAM SortedByCoordinate" ] ; then
		STAR_SAM_OUTPUT=$STAR_PREFIX"Aligned.sortedByCoord.out.bam"
	else 
		STAR_SAM_OUTPUT=$STAR_PREFIX"Aligned.out.sam"
	fi

	IS_GZIPPED_INPUT=`expr match "$R1_SUFFIX" '.*\(gz$\)'`
	# Specify read files command - default to gunzip, but if ungzipped input switch off
	if [ "$IS_GZIPPED_INPUT" == '' ]; then
		read_files_command=""
	fi
	if [ $paired -eq 1 ]; then
		echo "Mapping files: $R1_FILE and $R2_FILE"
	else 
		echo "Mapping file: $R1_FILE (single mode)"
	fi
	echo "Results will be stored in $OUT_DIR"
	#map using STAR:
	if [ ! -e "$STAR_SAM_OUTPUT" ]; 
	then
		STAR --runThreadN $threads --genomeDir $GENOME_DIR --genomeLoad $genome_load_option --readFilesIn $R1_FILE $R2_FILE $read_files_command --outFileNamePrefix $STAR_PREFIX $quant_mode $wiggle_output_command --outSAMmode $out_sam_mode --outSAMtype $out_sam_type --outSAMstrandField $out_sam_strand_field $out_reads_unmapped --limitBAMsortRAM $limit_bam_sort_ram $additional_star_options 
	fi
done
# remove shared genome from the memory - cleaning after mapping
STAR --genomeDir $GENOME_DIR --genomeLoad remove 
echo "Finished processing"
