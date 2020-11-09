#
# Download original Sarbecovirus sequences matching the alignment of Boni et al., and prepare associated 
# metadata for prediction script
#

import pandas
from Bio import SeqIO
from Bio import Entrez

SOURCE_ALIGMENT = "ExternalData/sarbecovirus/boni_et_al_NRR1_alignment.fas"
GB_OUTPUT_PATH = "ExternalData/sarbecovirus/sarbecovirus_raw.gb"
CSV_OUTPUT_PATH = "ExternalData/sarbecovirus/sarbecovirus_metadata.csv"


# Load aligment and extract accession numbers
alignment = SeqIO.parse(SOURCE_ALIGMENT, "fasta")

alignment_names = (seq.id for seq in alignment)
alignment_names = (name.replace("||", "|") for name in alignment_names)  # At least one name has a spurious double delimiter
accessions = {name.split("|")[3]: name for name in alignment_names}


# Download original genbank data
#  TODO: Currently dropping all GISAID sequences
ncbi_accessions = [acc for acc in accessions.keys() if not acc.startswith("EPI_")]

with Entrez.efetch(db="nucleotide", id=ncbi_accessions, rettype="gb", retmode="text") as fetch_handle:
    gb_seqs = fetch_handle.read()
    
with open(GB_OUTPUT_PATH, "wt") as out_file:
    out_file.write(gb_seqs)


# Create a matching metadata file
gb_seqs = SeqIO.parse(GB_OUTPUT_PATH, "gb")
seq_ids = {seq.name: seq.id for seq in gb_seqs}  # Above accession numbers do not include version, seq.id has format "accession.version"

metadata = []

for acc in ncbi_accessions:
    metadata.append({"Name": accessions[acc],
                     "SequenceID": seq_ids[acc]})

pandas.DataFrame(metadata).to_csv(CSV_OUTPUT_PATH, index=False)
