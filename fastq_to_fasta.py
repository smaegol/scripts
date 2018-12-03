#!/usr/bin/env python

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

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
    description='Convert fastq to fasta')

parser.add_argument('--input', dest='inputfile', action='store',
                    help='Input fastq file (required)', required=True)
parser.add_argument('--output', dest='outputfile', action='store',
                    help='Output fasta file (required)', required=True)
args = parser.parse_args()


with open(args.inputfile, "rU") as input_handle:
    with open(args.outputfile,"w") as output_handle:
        sequences = SeqIO.parse(input_handle, "fastq")
        count = SeqIO.write(sequences, output_handle, "fasta")

print("Converted %i records" % count)

