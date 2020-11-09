#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 15:23:39 2018

@author: nardus

- Download the genbank files matching accessions in '../InternalData/AllInternalData_Checked.csv'
- Merge these into a single fasta file
- Files are named using only the accession number (because some virus names
  in the dataset contain slashes)
"""

import os
import pandas as pd
from Bio import Entrez, SeqIO
from urllib.error import HTTPError

Entrez.tool = "ZoonosisPredictor pipeline - DownloadSequences.py"

Entrez.email = input("Enter an email address to use NCBI e-utils: ")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Get and check accessions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
InternalData = pd.read_csv("./InternalData/Final_Accessions_Unique_Spp.csv")

# Split accessions of segmented viruses:
Accessions = []

for row in InternalData.iterrows():
    if not row[1].isnull()['Genbank.accession']:  # Some entries do not have accession numbers (no sequences available)
        acc = row[1]["Genbank.accession"]
        Accessions += acc.split("; ")

# Remove trailing spaces (if any):
Accessions = [acc.rstrip() for acc in Accessions]

    
# Check for duplicates (idicative of a copy/paste error during data entry)
UniqueAccessions = set(Accessions)

if len(UniqueAccessions) != len(Accessions):
    Seen = []
    Duplicates = []
    
    for acc in Accessions:
        if acc in Seen:
            Duplicates.append(acc)
        else:
            Seen.append(acc)
    
    raise ValueError("Duplicate accession number(s) found: {}".format(Duplicates))




# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Download genbank files
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
os.makedirs("./ExternalData/Sequences", exist_ok=True)
os.chdir("./ExternalData/Sequences")

for ID in Accessions:
    FileName = "{}.gb".format(ID)
    
    if not os.path.isfile(FileName):
        try:
            QueryHandle = Entrez.efetch(db="nucleotide", id=ID, 
                                        rettype="gb", retmode="text")
        except HTTPError as Error:
            if Error.code == 400:  # Bad request
                raise ValueError("Accession number not found: {}".format(Accessions(ID)))
            else:
                raise
                
        with open(FileName, 'w') as WriteHandle:
            WriteHandle.write(QueryHandle.read())
        
    else:
        print("{} exists, not updated".format(FileName))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Combine to fasta
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Sequences = []

for ID in Accessions:
    FileName = "{}.gb".format(ID)
    SeqRec = SeqIO.read(FileName, format = "genbank")
    # Make sure sequences are named by accession only
    SeqRec.id = SeqRec.name = SeqRec.description = ID
    Sequences.append(SeqRec)
    

with open("CombinedSequences.fasta", 'w') as OutFile:
    SeqIO.write(Sequences, OutFile, "fasta")




