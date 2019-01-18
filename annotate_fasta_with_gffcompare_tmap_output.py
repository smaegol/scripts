#!/usr/bin/env python

#######################################################################################
###                                                                                 ###
###     Copyright (C) 2019  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
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

# script annotates input fasta file with the output fron Gffcompare (tmap file)
# adds reg_gene_id, ref_id and class code to the sequence id

from Bio import SeqIO
import pandas as pd
import argparse
from tqdm import *


parser = argparse.ArgumentParser(
    description='Convert fastq to fasta')
	
parser.add_argument('--fasta', dest='fasta', action='store',
                    help='Input fasta file (required)', required=True)
parser.add_argument('--tmap', dest='tmap_file', action='store',
                    help='Input tmap file from gffcompare (required)', required=True)
parser.add_argument('--output', dest='outputfile', action='store',
                    help='Output fasta file (required)', required=True)
args = parser.parse_args()

#read in annotations from gffcompare
tmap_annotations = pd.read_csv(args.tmap_file, sep="\t")
tmap_annotations.set_index("qry_id",inplace=True)

seq_dict = SeqIO.index(args.fasta,"fasta")
number_of_records=len(seq_dict)

#open input files:
with open(args.fasta,"r") as input_fasta, open(args.outputfile,"w") as output_file:
    for record in tqdm(SeqIO.parse(input_fasta,"fasta"),total=number_of_records,desc="Annotating fasta file",dynamic_ncols=True,position=1):
        annotation = tmap_annotations.loc[record.id,["ref_gene_id","ref_id","class_code"]] # get annotation for given sequence id
        if len(annotation.shape)==2: # if multiple hits per query, take only first one
            annotation = annotation.iloc[0]
        else:
            annotation = annotation
        new_sequence_id=annotation["ref_id"] + "|" + annotation["ref_gene_id"] + "|" + annotation["class_code"] + "|" + record.id # create new sequence id
        record.id=new_sequence_id

        SeqIO.write(record, output_file, "fasta")  # write sequence record to output
print("\n")
