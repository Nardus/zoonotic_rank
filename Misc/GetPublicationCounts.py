#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 14:16:14 2018

@author: nardus

- Get the number of publications matching each virus name in PubMed
- Search performed *at the species level* for all viruses listed in 
  InternalData/CuratedTaxonomy_MSL2018v1.csv
"""

import json
import csv
import time
import warnings
import pandas as pd
from Bio import Entrez

Entrez.tool = "ZoonosisPredictor pipeline - GetPublicationCounts.py"

Entrez.email = input("Enter an email address to use NCBI e-utils: ")
api_key = input("If you have an NCBI API key, enter it here (or leave empty): ")

if len(api_key) == 0:
    warnings.warn("No API key: Requests will be limited to 2 searches per second.")
    api_key = None


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Utils
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
class Searcher(object):
    """
     A basic implementation of esearch with automatic rate-limiting
    """
    
    def __init__(self, db, apiKey = None):
        self.requestsMade = 0
        self.lastRequest = 0
        
        self.db = db
        self.apiKey = apiKey
        
        if apiKey is not None:
            self.requestsPerSecond = 10
        else:
            self.requestsPerSecond = 2  # Since December 2018, only 2 requests/second
                                        # is allowed by NCBI (unless an API key is given)
                                        # (their docs say 3, but that doesn't seem to work...)
                                        
                                   
    def search(self, term):
        # check if we need to rate limit ourselves
        if self.requestsMade >= self.requestsPerSecond:
            elapsed = time.time() - self.lastRequest
            
            if elapsed < 1:
                time.sleep(1 - elapsed)
                
            self.requestsMade = 0
        
        # Send the query
        self.requestsMade += 1
        self.lastRequest = time.time()
        handle = Entrez.esearch(db=self.db, term=term, api_key = self.apiKey)
        record = Entrez.read(handle)
        
        return record
        



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Get virus names
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
virusData = pd.read_csv('./CalculatedData/FinalData_Cleaned.csv')

virusNames = set(virusData.LatestSppName)
alternativeNames = {row.LatestSppName: row.UniversalName for row in virusData.itertuples()}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Get publication counts
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
publicationCounts = {}

taxonomy = Searcher('taxonomy', api_key)
pubmed = Searcher('pmc', api_key)  # Only PubMedCentral can be searched by taxID...


for virus in virusNames:
    # Look up taxid:
    taxResult = taxonomy.search(virus)
    
    if taxResult['Count'] == '0':
        # Try an alternate name
        universalName = alternativeNames[virus]
        
        if virus != universalName:
            warnings.warn('Taxonomy search failed for {0}, trying {1}'.format(virus, universalName), Warning)
            taxResult = taxonomy.search(universalName)
    
    if taxResult['Count'] != '1':
        raise KeyError('Taxonomy search failed: {0} returned {1} results'.format(virus, taxResult['Count']))
    
    if len(taxResult['IdList']) != 1:
        raise KeyError('Taxonomy search failed: {0} returned {1} TaxIDs'.format(virus, taxResult['IdList']))
    
    taxID = taxResult['IdList'][0]
    
    #Search pubmed
    searchTerm = 'txid{0}[Organism:exp]'.format(taxID)  # This is the query generated when clicking through from NCBI taxonomy (in the subtree links column)
    result = pubmed.search(searchTerm)
    
    # Don't accept partial matches (all parts of name must be found):
    if 'ErrorList' in result.keys():
        if len(result['ErrorList']['PhraseNotFound']) > 0:
            result['Count'] = 0
    
    publicationCounts[virus] = int(result["Count"])



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Save result
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
outName = "./InternalData/PublicationCounts.csv"

with open(outName, 'w') as outFile:
    writer = csv.writer(outFile)
    writer.writerow(['VirusName', 'PubmedResults'])
    writer.writerows(publicationCounts.items())






