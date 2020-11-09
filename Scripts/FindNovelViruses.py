#!/usr/bin/env python3

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Detect recently added viruses in ICTV taxonomy and download sequences for them
#
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import argparse
import pandas
import warnings
import re
from Bio import Entrez
from Bio import SeqIO
from math import isnan

class SearchError(Exception):
    pass

parser = argparse.ArgumentParser(description="Detect recently added viruses in ICTV taxonomy and download sequences for them.")
parser.add_argument("msl", type=str, help="Path to an ICTV master species list")
parser.add_argument("vmr", type=str, help="Path to a matching virus metadata resource file")
parser.add_argument("output_base", type=str, help="Basename for output files. Output generated are {output_base}.gb (the sequences) and {output_base}.csv (matching metadata).")

inputargs = parser.parse_args()

Entrez.email = input("Enter an email address to use NCBI e-utils: ")
api_key = input("If you have an NCBI API key, enter it here to speed up downloads (or leave empty): ")

if len(api_key) == 0:
    warnings.warn("No API key: Number of requests will rate-limited")
    api_key = None

Entrez.api_key = api_key


# Families to consider:
# (these are listed as infecting animals in either Field’s Virology or the ViralZone database)
# NOTE: Some of these are not present in the training data: Anelloviridae, Genomoviridae
FAMILIES = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Arteriviridae", 
            "Asfarviridae", "Astroviridae", "Birnaviridae", "Bornaviridae",
            "Caliciviridae", "Circoviridae", "Coronaviridae", "Filoviridae",
            "Flaviviridae", "Genomoviridae", "Hantaviridae", "Hepadnaviridae",
            "Hepeviridae", "Herpesviridae", "Matonaviridae", "Nairoviridae",
            "Orthomyxoviridae", "Papillomaviridae", "Paramyxoviridae",
            "Parvoviridae", "Peribunyaviridae", "Phasmaviridae", "Phenuiviridae",
            "Picobirnaviridae", "Picornaviridae", "Pneumoviridae",
            "Polyomaviridae", "Poxviridae", "Reoviridae", "Retroviridae",
            "Rhabdoviridae", "Sunviridae", "Tobaniviridae", "Togaviridae"]

# Data-entry errors in 2019 VMR:
accession_replacements = {"741759": "KF741759",
                          "NC_129128": "NC_029128",
                          "EU7257772": "EU725772"}


# Find novel viruses
msl_raw = pandas.read_excel(inputargs.msl, sheet_name=2) # Read sheet number 3 (first 2 are metadata)
vmr_raw = pandas.read_excel(inputargs.vmr, sheet_name=0)

msl_raw.columns = [c.replace(' ', '_') for c in msl_raw.columns]
vmr_raw.columns = [c.replace(' ', '_') for c in vmr_raw.columns]

msl_subset = msl_raw[msl_raw.Family.isin(FAMILIES)]
novel_viruses = msl_subset[msl_subset.Last_Change.str.contains("New,")]

# Keep only those for which a complete genome is available
accession_data = vmr_raw[vmr_raw.Species.isin(novel_viruses.Species)]
accession_data = accession_data[accession_data.Genome_coverage.notna()]
accession_data = accession_data[accession_data.Genome_coverage.str.match("Complete genome")]

# Parse accessions
accessions = {}

for row in accession_data.itertuples():
    if not isinstance(row.Virus_REFSEQ_accession, str):
        assert isnan(row.Virus_REFSEQ_accession)
        acc_string = row.Virus_GENBANK_accession
    else:
        acc_string = row.Virus_REFSEQ_accession
    
    # Split accessions of segmented viruses:
    acc_list = acc_string.split(";")
    acc_list = [a.lstrip() for a in acc_list]
    
    if len(acc_list) > 1:
        acc_list = [re.search("[\S]+:[\s]*([\S]+)", a).group(1) for a in acc_list]
    
    for acc in acc_list:
        # Fix data-entry errors in 2019 VMR:
        if acc in accession_replacements.keys():
            warnings.warn("Invalid accession '{0}' replaced with '{1}'".format(acc, accession_replacements[acc]))
            acc = accession_replacements[acc]
        
        accessions[acc] = row.Species


# Download sequences
with Entrez.efetch(db="nucleotide", id=list(accessions.keys()), rettype="gb", retmode="text") as fetch_handle:
    gb_seqs = fetch_handle.read()

seq_output_path = "{}.gb".format(inputargs.output_base)

with open(seq_output_path, "wt") as out_file:
    out_file.write(gb_seqs)


# Create a matching metadata file
gb_seqs = SeqIO.parse(seq_output_path, "gb")
seq_ids = {seq.name: seq.id for seq in gb_seqs}  # Above accession numbers do not include version, seq.id has format "accession.version"

metadata = []

for acc in accessions.keys():
    metadata.append({"Name": accessions[acc],
                     "SequenceID": seq_ids[acc]})

csv_output_path = "{}.csv".format(inputargs.output_base)
pandas.DataFrame(metadata).to_csv(csv_output_path, index=False)