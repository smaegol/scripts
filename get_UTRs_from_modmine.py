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
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from intermine.webservice import Service
import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description="Get worm' UTRs sequences from Modmine")

parser.add_argument('--utr', dest='utr', action='store',
                    help='UTR to collect, default = 3utr', choices = ['5utr','3utr'], default="3utr", required=False)
parser.add_argument('--output', dest='output_fasta', action='store',
                    help='Output fasta file (required)', required=True)
parser.add_argument('--output_tsv', dest='output_tsv', action='store',
                    help='Output tsv file with results (optional)', required=False)
parser.add_argument('--gene', dest='gene_symbol', action='store',
                    help='Gene we want UTR sequence for (symbol) (optional)', required=False)
parser.add_argument('--source', dest='annotation_source', action='store',
                    help='Source of annotations (optional) default = WormBase', choices = ["WormBase","CEUP1", "UTRome_V2"], default="WormBase", required=False)
parser.add_argument('--no_hits', dest='no_hits', action='store',
                    help='Limit number of hits to specified number',required=False)
parser.add_argument('--utrs_only', dest='utrs_only', action='store_true',
                    help='Return only genes which have UTR',required=False)
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

# available sources dict
available_sources = {"WormBase": "WormBase", "CEUP1" : "modENCODE-3'UTRome CEUP1 sequences, alignments, annotation","UTRome_V2" : "modENCODE-UTRome_V2_3UTRs_multiple_evidences"}

#connect to ModMine service
service = Service("http://intermine.modencode.org/release-33/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")
# Type constraints should come early - before all mentions of the paths they constrain

# Specify which UTRs to collect
if (args.utr=='3utr'):
    query.add_constraint("UTRs", "ThreePrimeUTR")
elif (args.utr=='5utr'):
    query.add_constraint("UTRs", "FivePrimeUTR")
else:
    print("Wrong UTR option provided. Exiting...")
    exit()

# if specified, limit results for given gene
if (args.gene_symbol):
    query.add_constraint("symbol", "=", args.gene_symbol)


# The view specifies the output columns
query.add_view(
    "source", "UTRs.primaryIdentifier", "UTRs.length",
    "UTRs.chromosome.primaryIdentifier", "UTRs.chromosomeLocation.start",
    "UTRs.chromosomeLocation.end", "primaryIdentifier", "sequence.length",
    "UTRs.sequence.residues"
)

# Select annotation source
if (str(args.annotation_source) in available_sources.keys()):
    annotation_source = available_sources[args.annotation_source]
else:
    print ("Wrong annotation source provided.Exiting...")
    exit()
query.add_constraint("source", "=", annotation_source, code = "B")

# By default will show all genes, even if they don't have annotated UTRs
# When utrs_only option is on will not perform outerJoin and will show only genes with UTRs
if (args.utrs_only==False):
    query.outerjoin("UTRs")

#Get number of hits returned by query
number_of_hits = query.count()

if (number_of_hits > 0):
    print ("Found " + str(number_of_hits) + " hits in the database")
else:
    print ("no results found. Try to change your input options")
    exit()

# Introduce limit on hits (if specified)
if (args.no_hits):
    limit_hits = args.no_hits
else:
    # if not specified, use all hits
    limit_hits = number_of_hits

if (args.output_tsv):
    output_table = pd.DataFrame() # pandas data frame to store all results
with open(args.output_fasta,"w") as output_file:
    for row in query.results(size=10,row="dict"):
        if (args.output_tsv):
            results_pd = pd.Series(row)
            output_table = output_table.append(results_pd,ignore_index=True)
        seqrec = SeqRecord(Seq(str(row["Gene.UTRs.sequence.residues"])),row["Gene.primaryIdentifier"],description=str(row["Gene.symbol"]) + " " + str(row["Gene.UTRs.primaryIdentifier"]))
        SeqIO.write(seqrec, output_file, "fasta")

print ("Finished writing sequences to " + args.output_fasta)
if (args.output_tsv):
    print("Writing summary table to " + args.output_tsv + "...")
    output_table.to_csv(args.output_tsv,sep="\t")
