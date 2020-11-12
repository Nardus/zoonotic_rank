#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract gene sequences and calculate genomic features by calling Richard Orton's
CPB_Machine.jar


@author: nardus
"""

import argparse
import pandas
import tempfile
import os
import shutil
import warnings
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from subprocess import run


# --------------------------------------------------------------------------------------------------------------------------
# Resolve relative location of feature calculation tools
# - CPB_Machine.jar expected to be in a subdirectory called 'External' relative to this script
# - Script for non-coding sequences expected in the same directory as this script
# --------------------------------------------------------------------------------------------------------------------------
scriptPath = os.path.dirname(os.path.realpath(__file__))
CPB_LOCATION = os.path.join(scriptPath, 'External', 'CPB_Machine.jar')

NC_LOCATION = os.path.join(scriptPath, 'calc_noncoding_genome_features.py')


# --------------------------------------------------------------------------------------------------------------------------
# Get the current path so we can return output here
# --------------------------------------------------------------------------------------------------------------------------
START_PATH = os.getcwd()



# --------------------------------------------------------------------------------------------------------------------------
## Parse input:
# --------------------------------------------------------------------------------------------------------------------------
argParser = argparse.ArgumentParser(description='Extract gene sequences and calculate genomic features')

argParser.add_argument('metadata', metavar = 'metadata', type = str, 
                       help = 'Path to a meta-data file with columns "Name" and "SequenceID", in csv format. Names can be repeated for different SequenceIDs, in which case a single summary will be produced.')

argParser.add_argument('sequences', metavar = 'sequences', type = str,
                       help = 'Location of the sequence file(s) - can be a single file or a directory containing multiple files. If option "--extract" is used, one or more GenBank flat file(s); otherwise, fasta file(s). If a directory is supplied, files should have names in the format [SequenceID].gb or [SequenceID].fasta')

argParser.add_argument('--extract', action = 'store_true',
                       help = 'Extract coding sequences from individual GenBank files. When not used, sequences are expected to be in fasta format.')

argParser.add_argument('--noncoding', action = 'store_true',
                       help = 'Treat all sequences as non-coding. In this case only a subset of features will be returned (since codon bias, etc. do not apply).')

argParser.add_argument('--out', metavar='output_file', type=str,
                       help = 'Output file location. If not specified, results will be named [metadata]_GenomeFeatures.csv, where "metadata" is the full path to the metadata file speciefied above.')

inputArgs = argParser.parse_args()


## Check input
# Check input file / folder exists:
if not os.path.exists(inputArgs.sequences):
    raise FileNotFoundError('{} not found'.format(inputArgs.sequences))



# --------------------------------------------------------------------------------------------------------------------------
## Functions:
# --------------------------------------------------------------------------------------------------------------------------
def fuzzy_extract(dictionary, key, returnType = 'string'):
    """Extract a value from a dictionary which may not contain the key.
    Instead of failing when key is not present, returns either an empty
    string or 0, depending on type requested
    """
    
    try:
        result = dictionary[key]
    except KeyError:
        if returnType == 'string':
            result = ''
        elif returnType == 'intList':
            result = [0]
        else:
            raise NotImplementedError('Invalid returnType')
            
    return result


def check_cds(feature, extractedSeq, name = None, problems_file = None):
    """Check extracted CDS sequences for validity, returning an informative warning 
    messages before skipping/discarding invalid ones.
    Data on skipped sequences will be appended to the open file handle 'problems_file',
    if supplied (in csv format with columns: name, accession, protein_name, protein_id, length, problem_description).
    """
    
    accession = extractedSeq.id
    proteinName = ", ".join(fuzzy_extract(feature.qualifiers, "product"))
    proteinID = ", ".join(fuzzy_extract(feature.qualifiers, "protein_id"))
    baseMessage = "Sequence {0}: {1} (protein id {2})".format(accession, proteinName, proteinID)
    
    # Check if for causes of concern:
    wrongLength = len(extractedSeq) % 3 != 0
    truncatedStart = type(feature.location.start) == BeforePosition
    truncatedEnd = type(feature.location.end) == AfterPosition
    
    startPos = fuzzy_extract(feature.qualifiers, "codon_start", "intList")[0]
    shiftedStart = int(startPos) > 1
    
    
    # Respond to these    
    if truncatedStart or truncatedEnd:
        problem = "truncation"
        message = "CDS skipped - truncation: {}".format(baseMessage)
        
    elif shiftedStart:
        problem = "translation start occurs before sequence start (truncation?)"
        message = "CDS skipped - shifted translation start found (codon {0}), assuming truncation: {1}".format(startPos, baseMessage)
        
    elif wrongLength:
        problem = "length not divisible by 3 (cause not identified)"
        message = "CDS skipped - length not a multiple of 3, but cause not identified: {}".format(baseMessage)
        
    else:
        problem = ""
        message = ""
        
    
    # Return:
    if len(message) != 0:
        warnings.warn(message, Warning)
        
        if problems_file is not None:
            outline = ",".join([name, accession, proteinName, proteinID, str(len(extractedSeq)), problem])
            outline = "{0}\n".format(outline)
            _ = problems_file.write(outline)
        
        return None
    
    return extractedSeq


def rename_sequence(seqrecord, name):
    """Rename a seqrecord object so the specified name gets used when 
    writing to fasta
    """
    seqrecord.name = name
    seqrecord.id = name
    seqrecord.description = ''
    
    return seqrecord 


def extract_and_name_cds(seqrecord, name, problems_file = None):
    """Extract all (valid) CDS sequences from a parsed genbank record and rename
    them to a common name.
    If given, 'problems_file' should be an open file handle to which data on problematic 
    sequences can be recorded.
    """
    codingSeqs = []
    
    for feature in seqrecord.features:
        if feature.type.lower() == "cds":
            seq = feature.extract(seqrecord)
            validSeq = check_cds(feature, seq, name, problems_file)
            
            if validSeq is not None:
                renamedSeq = rename_sequence(validSeq, name)
                codingSeqs.append(renamedSeq)
    
    return codingSeqs
    



# --------------------------------------------------------------------------------------------------------------------------
## Read input data
# --------------------------------------------------------------------------------------------------------------------------
metaData = pandas.read_csv(inputArgs.metadata)

seqNames = {name: [] for name in set(metaData.Name)}

for index, rowData in metaData.iterrows():
    seqNames[rowData.Name].append(rowData.SequenceID)


# Load sequence data
if inputArgs.extract:
    fileFormat = "genbank"
    basename = "{}.gb"
else:
    fileFormat = "fasta"
    basename = "{}.fasta"


if os.path.isdir(inputArgs.sequences):
    # Load data from individual files in this folder:
    os.chdir(inputArgs.sequences)
    seqData = {}
        
    for name in seqNames.keys():
        for acc in seqNames[name]:
            filename = basename.format(acc)
            seqData[acc] = SeqIO.read(filename, format = fileFormat)
else:
    # A single file:
    seqData = SeqIO.parse(inputArgs.sequences, format = fileFormat)
    seqData = {seq.id: seq for seq in seqData}




# --------------------------------------------------------------------------------------------------------------------------
## Create a fasta file with gene sequences:
#  - If needed, all CDS sequences extracted from each GenBank file
#  - Sequences renamed to have a common name (so CPB_Machine.jar) creates a single summary for each 'Name' in the meta data
#  - Immediately write to fasta
# --------------------------------------------------------------------------------------------------------------------------

## Create a temporary fasta file
#  - To write a single-line fasta, we have to use SeqIO.FastWriter directly:
outFasta = tempfile.NamedTemporaryFile("wt", delete = False, suffix = ".fasta")

outWriter = SeqIO.FastaIO.FastaWriter(outFasta, wrap=None)
outWriter.write_header()

#  - If we're extracting CDS regions, record data on problematic coding regions too:
if not inputArgs.noncoding:
    problemFile = tempfile.NamedTemporaryFile("wt", delete = False, suffix = ".csv")
    problemFile.write("name,accession,protein_name,protein_id,length,problem_description\n")
else:
    problemFile = None


# Process the sequence data:
processedIDs = []

for name in seqNames.keys():
    newName = "....{}_cds_".format(name)
        
    for seqID in seqNames[name]:
        if inputArgs.extract and not inputArgs.noncoding:
            sequences = extract_and_name_cds(seqData[seqID], newName, problems_file = problemFile)
        else:
            sequences = [rename_sequence(seqData[seqID], newName)] # Expect a single seqrecord
            
            # Do a basic check to ensure sequence is valid:
            if not inputArgs.noncoding and len(sequences[0]) % 3 != 0:
                warnings.warn("Sequence skipped ({}): Length is not divisble by 3. Use the --noncoding option to calculate features for noncoding sequences.".format(name), Warning)
                problem = "length not divisible by 3 (cause not identified)"
                outline = ",".join([newName, "", "", "", str(len(sequences[0])), problem])
                outline = "{0}\n".format(outline)
                _ = problemFile.write(outline)
                sequences = None
        
        if sequences is not None:
            outWriter.write_records(sequences)
        
        processedIDs.append(seqID)
        

outWriter.write_footer()
outFasta.close()

if problemFile is not None:
    problemFile.close()


if not inputArgs.extract:
    notProcessed = [key for key in seqData.keys() if key not in processedIDs]
    
    if len(notProcessed) != 0:
        message = "Fasta file contained {} sequences that were not present in the meta-data".format(len(notProcessed))
        warnings.warn(message, Warning)



## Create an 'accessions file'
#  - CPB_Machine does not use the last column (except for parsing)
outData = tempfile.NamedTemporaryFile("wt", delete = False, suffix = ".csv")

outDF = pandas.DataFrame()
outDF["Name"] = seqNames.keys()
outDF["TaxID"] = "NA"
outDF["Species"] = "NA"
outDF["Family"] = "NA"

outDF.to_csv(outData, sep='\t', index = False, header = False)
outData.close()


# --------------------------------------------------------------------------------------------------------------------------
# Calculate the features
# --------------------------------------------------------------------------------------------------------------------------
outputPath = os.path.splitext(outFasta.name)[0]
outputFile = "{}_dat.txt".format(outputPath)  # This path is used by CPB_Machine by default


if not inputArgs.noncoding:
    # Use CPB_Machine if possible
    run(["java", "-jar", CPB_LOCATION, outFasta.name, outData.name], check=True)
    
else:
    # In this case the sequence lengths might not be divisible by 3, so use 
    #  the python version to calculate relevant features only
    run(["python3", NC_LOCATION, outFasta.name, outputFile], check=True)



# --------------------------------------------------------------------------------------------------------------------------
# Write final output and clean up temp files
# --------------------------------------------------------------------------------------------------------------------------
inputPath = os.path.splitext(inputArgs.metadata)[0]

if inputArgs.out is None:
    # Copy results to same directory as input data
    finalOutName = "{}_GenomeFeatures.csv".format(inputPath)
    finalProbsName = "{}_Problems.csv".format(inputPath)
else:
    outbase = os.path.splitext(inputPath)[0]  # Full path, minus the extension
    finalOutName = inputArgs.out
    finalProbsName = "{}_Problems.csv".format(outbase)

os.chdir(START_PATH)
shutil.copy(outputFile, finalOutName)  # Need to copy then remove, in case files are on different partitions (in which case os.rename won't work)
os.remove(outputFile)

if problemFile is not None:
    shutil.copy(problemFile.name, finalProbsName)
    os.remove(problemFile.name)

# Remove temporary files
os.unlink(outFasta.name)
os.unlink(outData.name)
