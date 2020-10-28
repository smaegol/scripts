#!/bin/bash

## Processing direct RNA seq by nanopore with Nanopolish

#######################################################################################
###                                                                                 ###
###     Copyright (C) 2020  Pawel Krawczyk (pkrawczyk@iimcb.gov.pl)                 ###
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




# define parameters
threads=20
nanopolish_output_folder="nanopolish"
guppy_folder="guppy_out_3.6.0"
merged_reads='merged_reads.fastq'
nanopolish_index=$merged_reads".index"
seqkit_chunks=1
parallel_nanopolish=4
genome='human'

#load all required modules
module load minimap2
module load samtools/1.9
module load nanopolish
module load bedtools

#Get variables from command line
while getopts ":n:t:g:s:o:r:h::" optname
  do
    case "$optname" in
        "h")
                echo "Nanopore analysis pipeline"
                echo ""
                echo "Options :"
                echo " -n : sample name (Required)"
                echo " -t: threads"
                echo " -o: nanopolish output folder"
                echo " -g: guppy output folder"
		echo " -s: split input into n chunks [provide as an argument to s] using seqkit"
		echo " -r: reference to use for Nanopolish [human|mouse|yeast|MHV|arabidopsis|worm|covid|pombe]"
                echo " -h : THIS HELP"
                echo ""
                exit 1
                ;;
        "n")
                sample_name=$OPTARG
                ;;
	"s")
                seqkit_chunks=$OPTARG
                ;;
        "t")
                threads=$OPTARG
                ;;
        "o")
                nanopolish_output_folder=$OPTARG
                ;;
        "g")
                guppy_folder=$OPTARG
                ;;
        "r")
                genome=$OPTARG
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

if [  -e $sample_name ]; then
	echo "Missing sample name [-n] (required)"
	exit 1
fi

# Define references
if [ $genome == 'human' ]; then
	reference_fasta="/home/smaegol/data/indexes/Hsapiens/annotation/gencode_32/gencode.v32.transcripts.with_ribosomal_from_full.fasta"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_hsapiens_gencode32_prompts_rna"
elif [ $genome == 'mouse' ]; then
	reference_fasta="/home/smaegol/data/indexes/Mus/annotation/gencode.vM22.transcripts.fa"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_mouse_gencode22"
elif [ $genome == 'yeast' ]; then
	reference_fasta="/home/smaegol/data/indexes/Scerevisiae/new_annot/proper_one/SacCer3_custom_annotation_ORFs_ncRNA_CUTs_SUTs_XUTs.fasta"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_yeast_new_PK_annot_CUTs_XUTs_ncRNA"
elif [ $genome == 'MHV' ]; then
	reference_fasta="/home/smaegol/data/indexes/coronaviruses/MHV/MHV.fasta"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_MHV"
elif [ $genome == 'arabidopsis' ]; then
	reference_fasta="/home/smaegol/data/indexes/Athaliana/TAIR10/annotation/Arabidopsis_thaliana.TAIR10.cdna.ncrna.all.fa"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_Athaliana"
elif [ $genome == 'worm' ]; then
	reference_fasta="/home/smaegol/data/indexes/Celegans/annotation/Caenorhabditis_elegans.WBcel235.cdna_ncrna.all.fa"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_celegan"
elif [ $genome == 'pombe' ]; then
	reference_fasta="/home/smaegol/data/indexes/SPombe/cds_introns_utrs.fa"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_spombe"
elif [ $genome == 'covid' ]; then
	reference_fasta="/home/smaegol/data/indexes/coronaviruses/SARS-Cov-2/MT007544.1.fasta"
	reference_minimap2=$reference_fasta".mmidx"
	analysis_suffix="_sars_cov_2"
else

	echo "Unknown reference type provided. Quiting..."
	exit 1
fi



current_folder=`pwd`

# Print input parameters

echo -e "\e[1mProvided parameters:\e[0m"
echo -e "sample_name: $sample_name"
echo -e "threads to use: $threads"
echo -e "output folder: $nanopolish_output_folder"
echo -e "guppy output folder: $guppy_folder"
echo -e "current folder: $current_folder"
echo -e "split reads into chunks: $seqkit_chunks"


bam_output=$sample_name'_'$analysis_suffix'.bam'
nanopolish_output=$sample_name'_'$analysis_suffix'_nanopolish.tsv'


if [ ! -d $nanopolish_output_folder ]; then 
	mkdir $nanopolish_output_folder
fi



# merge reads from guppy output

echo -e "\e[1mMerging reads\e[0m"
if [ -e $nanopolish_output_folder/$merged_reads ]; then
	echo -e "Reads seem to be already merged. Skipping..."
else
	cat $guppy_folder/*.fastq > $nanopolish_output_folder/$merged_reads
fi

#cd $nanopolish_output_folder

# make symbolic links of workspace and sequencing_summary.txt - required for nanopolish index and polya functions 
if [ -e $nanopolish_output_folder/workspace ]; then
	echo -e "Workspace is already present in the nanopolish folder"
else 
	ln -s $guppy_folder/workspace $nanopolish_output_folder/
fi

if [ -e $nanopolish_output_folder/sequencing_summary.txt ]; then
	echo -e "Sequencing_summary is already present in the nanopolish folder"
else 
	ln -s $guppy_folder/sequencing_summary.txt $nanopolish_output_folder/
fi
 
# Index reads with Nanopolish

cd $nanopolish_output_folder

echo -e "\e[1mIndexing reads with nanopolish\e[0m"
if [ ! -e $nanopolish_index ]; then
	nanopolish index -d workspace -s sequencing_summary.txt $merged_reads
else
	echo -e "Nanopolish index is already present. Skipping indexing..."
fi

if [ $seqkit_chunks -gt 1 ]; then
# if reads have to be splitted using seqkit
	module load seqkit
	echo -e "\e[1mSplitting reads into $seqkit_chunks chunks to improve processing time\e[0m"
	if [ ! -d $merged_reads".split" ]; then

		seqkit split -p $seqkit_chunks $merged_reads
	else
		number_fastq_files_present=`find merged_reads.fastq.split/ -name "*fastq" | wc -l`
		if [ $number_fastq_files_present == $seqkit_chunks ]; then
			echo -e "Sequences seem to be already splitted into chunks. Skipping seqkit split part"
		else 
			echo -e "Number of fastq files in $merged_reads.split folder does not equal specified number of chunks ($seqkit_chunks). Deleting and rerunning..."
			rm -rf $merged_reads".split"
			seqkit split -p $seqkit_chunks $merged_reads
		fi
	fi 
	
	echo -e "\e[1mparallel mapping with minimap2\e[0m"
	number_bam_files_present=`find merged_reads.fastq.split/ -name "*$analysis_suffix.bam" | wc -l`
	if [ $number_bam_files_present == $seqkit_chunks ]; then
		echo -e "$number_bam_files_present bam files already present in the folder. Skipping mapping..."
	else
		find . -name "*.fastq" | parallel -j 2 minimap2 -ax map-ont -t $threads --secondary=no $reference_minimap2 {} '|' samtools view -b -F 2320 -@ $threads '|' samtools sort -o {}"_$analysis_suffix.bam" -@ $threads
	fi

	find . -name "*$analysis_suffix.bam" | parallel -j $threads samtools index {}


	bedtools bamtobed -cigar -i $merged_reads"_"$analysis_suffix".bam" > "merged_reads.fastq__"$analysis_suffix".bam.bed"

	echo -e "\e[1mestimating poly(A) length in parallel using $parallel_nanopolish simultaneous Nanopolish threads\e[0m"
	
	find $merged_reads".split"/ -name "*$analysis_suffix.bam" | parallel -j $parallel_nanopolish nanopolish polya -t $threads -b {} -r $merged_reads -g $reference_fasta '>' {}'_'$analysis_suffix'_nanopolish.tsv' 

	echo -e "\e[1mFinished parallel Nanopolish processing. Merging predictions\e[0m".

	cat */*$analysis_suffix'_nanopolish.tsv' > "tmp."$sample_name"_"$analysis_suffix"_nanopolish.tsv"
	head -n 1 "tmp."$sample_name"_"$analysis_suffix"_nanopolish.tsv" > $sample_name"_"$analysis_suffix"_nanopolish.tsv"
	tail -n +2 "tmp."$sample_name"_"$analysis_suffix"_nanopolish.tsv" | grep -v "^read" >> $sample_name"_"$analysis_suffix"_nanopolish.tsv"

	echo  -e "\e[1mFinished...\e[0m"
else
	echo -e "\e[1mProcessing single file - no splitting into chunks\e[0m"
	echo -e "\e[1mmapping reads with minimap2\e[0m"
	if [ ! -e $bam_output ]; then
		minimap2 -ax map-ont -t $threads --secondary=no $reference_minimap2 $merged_reads | samtools view -b -@ $threads -F 2320 | samtools sort -@ $threads -o $bam_output
		echo -e "\e[1mindexing bam file with samtools\e[0m"
		samtools index $bam_output
		echo -e "\e[1mProducing bed output\e[0m"
		bedtools bamtobed -cigar -i $bam_output > $bam_output".bed"
	else
		echo -e "Bam file with mapping is already present. Skipping mapping..."
	fi
	echo -e "\e[1mEstimating polyA length\e[0m"
	nanopolish polya -t $threads -g $reference_fasta -r $merged_reads -b $bam_output > $nanopolish_output
fi
echo -e "\e[1mFinished\e[0m"
