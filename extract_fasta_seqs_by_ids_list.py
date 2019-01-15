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

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='Convert fastq to fasta')
	
parser.add_argument('--fasta', dest='inputfasta', action='store',
                    help='Input fasta file (required)', required=True)
parser.add_argument('--ids', dest='ids_list', action='store',
                    help='Ids list (required)', required=True)
parser.add_argument('--output', dest='outputfile', action='store',
                    help='Output fasta file (required)', required=True)
args = parser.parse_args()

print("Reading sequence ids...")

# Read sequence ids from the file
with open(args.ids_list) as ids_file:
    seq_ids = ids_file.read().splitlines() 

print("Reading and filtering sequences ...")
    
# Read fasta file and output filtered sequences
with open(args.inputfasta,"r") as input_fasta, open(args.outputfile,"w") as output_file:
    for record in SeqIO.parse(input_fasta,"fasta"):
        if record.id in seq_ids:
            SeqIO.write(record, output_file, "fasta")

print("Finished")