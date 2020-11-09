#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 15:23:39 2018

@author: nardus

- Download genbank files matching ENSEMBL IDs in each of the human gene set files specified below
- Files are named using only the ENSEMBL ID
"""

import requests
import time
import json
import warnings
import os
import pandas
from time import sleep


DATA_PATH = os.path.join('..', 'CalculatedData', 'HumanGeneSets')
ISG = os.path.join('..', 'InternalData', 'Shaw2017_raw', 'ISG_CountsPerMillion_Human.csv')
HOUSEKEEPING = os.path.join('..', 'ExternalData', 'HousekeepingGenes.txt')


# --------------------------------------------------------------------------------------------------------------------------
## ENSEMBL API
# -------------------------------------------------------------------------------------------------------------------------- 
class Ensembl(object):
    """A class for managing requests to Ensembl's REST API"""
    
    def __init__(self, server = 'https://rest.ensembl.org', requestsPerSecond = 15, timeout = 10000, chunksize = 50):
        self.server = server
        self.requestsPerSecond = requestsPerSecond
        self.timeout = timeout  # only used for post requests, where reading response may take a while
        self.requestsMade = 0
        self.lastRequest = 0
        
        if chunksize > 1000:  # this is the maximum allowed by Ensembl (I think, though in reality 50 seems to be the max)
            warnings.warn('Chunksize set too high - falling back to 1000')
            self.chunksize = 1000
        else:
            self.chunksize = chunksize
    
    
    
    # Generic requester:
    def _ensembl_request(self, accessPoint, identifier, filters, idName = 'ids', post = False, contentType = 'application/json'):
        
        """Send a request to ENSEMBL and return raw response text
        If post = True, multiple identifiers are accepted. 
            - This will be more efficient than individual requests (which are 
              rate-limited)
            - But not all access points accept POST requests
        """
        
        # Build query
        if post:
            postString = '{0}/{1}'.format(self.server, accessPoint)
            filters[idName] = identifier
            postData = json.dumps(filters)
            
        else:
            filterStrings = ['{0}={1}'.format(key, value) for key, value in filters.items()]
            filterStrings  = ';'.join(filterStrings)
            getString = '{0}/{1}/{2}?{3}'.format(self.server, accessPoint, identifier, filterStrings)
        
        
        # check if we need to rate limit ourselves
        if self.requestsMade >= self.requestsPerSecond:
            elapsed = time.time() - self.lastRequest
            
            if elapsed < 1:
                time.sleep(1 - elapsed)
                
            self.lastRequest = time.time()
            self.requestsMade = 0
        
        # Send the query
        self.requestsMade += 1
        
        if post:
            response = requests.post(postString, data = postData,
                                     headers = {'Content-Type': contentType},
                                     timeout = self.timeout)  # Reading response may take a long time
            
        else:
            response = requests.get(getString, 
                                    headers = {'Content-Type': contentType},
                                    timeout = self.timeout)
        
        # Check response and return
        if not response.ok:
            if response.status_code == 429: # We are being rate limited
                waitingTime = response.headers['Retry-After']
                time.sleep(float(waitingTime))
                
            elif response.status_code == 400: # ID not found (Bad Request)
                result = json.loads(response.text)
                if 'error' in result.keys():
                    warnings.warn(result['error'])
                    return None
                    
                else: # This shouldn't happen
                    response.raise_for_status()
                
            else:
                response.raise_for_status()
                
        return response.text
    
    
    # POST requests:
    def _chunk_generator(self, l, n):
        """Split a list into chunks of size n
        From https://stackoverflow.com/a/1751478 """
        n = max(1, n)  # For last item remaining if list is not a multiple of n
        return (l[i:i+n] for i in range(0, len(l), n))
    
    
    def _post_request(self, accessPoint, identifiers, filters, idName = 'ids'):
        """Construct and send a POST request
        Handles the 1000 item limit of ENSEMBL POST by splitting the query if
          needed (https://github.com/Ensembl/ensembl-rest/wiki/POST-Requests)
        Returns a list of raw text/json responses """
        
        idChunks = self._chunk_generator(identifiers, self.chunksize)
        responses = []
        
        for ids in idChunks:
            result = self._ensembl_request(accessPoint, ids, filters, 
                                           idName = idName,
                                           post = True)
            responses.append(result)
        
        return responses


    # Specific actions:
    def _get_transcripts(self, stableIDs, fromGeneID = True):
        """Get all transcripts associated with a list of ENSEMBL ID
        Returns a dictionary of dictionaries: {stableID: <..transcriptdata..>}"""
        
        accessPoint = "lookup/id"
        
        filters = {'species': 'human'}
        
        if fromGeneID:
            filters['expand'] = 1  # Ensures response will include nested info (e.g. transcripts)
        
        responses = self._post_request(accessPoint, stableIDs, filters)
        results = {}
        
        for resp in responses:
            resp = json.loads(resp)
            results.update(resp)
        
        return results
        
    
    def _check_ccds(self, transcriptID):
        """Check whether a transcript is in the CCDS database """
        crossrefs = self._ensembl_request('xrefs/id', transcriptID, {})
        crossrefs = json.loads(crossrefs)
        
        databases = [refdict['dbname'] for refdict in crossrefs]
        
        return ('CCDS' in databases)
    
    
    def _check_protein_coding(self, transcriptIDs):
        """Check whether transcripts are protein-coding
            Expects a list of transcriptIDs
            Returns a dictionary: {TranscriptID: <True/False>} """
        transcriptData = self._get_transcripts(transcriptIDs, fromGeneID = False)
        result = {}
        
        for key, value in transcriptData.items():
            if value is None:
                result[key] = None  # ID not found
            else:
                result[key] =  value['biotype'] == 'protein_coding'
        
        return result


    
    def _get_cds(self, transcriptIDs):
        """Get the coding sequences associated with transcript ids
            Returns a list of dictionaries """
            
        accessPoint = 'sequence/id'
        
        filters = {'type': 'cds',
                   'species': 'human'}
        
        responses = self._post_request(accessPoint, transcriptIDs, filters)
        results = []
        
        for resp in responses:
            results += json.loads(resp)
        
        return results
    
    
    
    def _get_cDNA(self, transcriptIDs):
        """Like _get_cds, but request a cDNA sequence instead (i.e. the coding
        sequence including any UTR's, but excluding introns, if present) """
        
        accessPoint = 'sequence/id'
        
        filters = {'type': 'cdna',
                   'species': 'human'}
        
        responses = self._post_request(accessPoint, transcriptIDs, filters)
        results = []
        
        for resp in responses:
            results += json.loads(resp)
        
        return results
 
        
    def summarise_transcripts(self, idList, fromGeneID = True):
        """Summarise the properties of transcripts 
            - When fromGeneID is true, all transcripts associated with the input 
              Ensembl gene IDs 
            - Otherwise, simply the properties associated with the input TranscriptIDs
        Returns a pandas dataframe """
        
        if type(idList) == str:
            idList = [idList]
        
        # Get transcripts and their properties
        transcriptData = self._get_transcripts(idList, fromGeneID)
        
        # Extract relevant data
        summary = {'GeneID': [],
                   'TranscriptID': [],
                   'Biotype': [],
                   'Canonical': [],
                   'InCCDS': []}
        
        for identifier in idList:
            if transcriptData[identifier] is None:
                # ID not found
                tData = None
            elif fromGeneID:
                tData = transcriptData[identifier]['Transcript']
            else:
                tData = [transcriptData[identifier]]
            
            # Add data (or indicate missing data):
            if tData is not None:
                for transcript in tData:
                    summary['GeneID'].append(transcript['Parent'])
                    summary['TranscriptID'].append(transcript['id'])
                    summary['Biotype'].append(transcript['biotype'])
                    summary['Canonical'].append(transcript['is_canonical'])
                    
                    inCCDS = self._check_ccds(transcript['id'])
                    summary['InCCDS'].append(inCCDS)
            else:
                if fromGeneID:
                    summary['GeneID'].append(identifier)
                    summary['TranscriptID'].append(None)
                else:
                    summary['TranscriptID'].append(identifier)
                    summary['GeneID'].append(None)
                
                summary['Biotype'].append(None)
                summary['Canonical'].append(None)
                summary['InCCDS'].append(None)
        
        return pandas.DataFrame(summary)
    
    
    def get_sequences(self, transcriptData, cds = True):
        """ Retrieve cds or cDNA sequences for all genes and transcripts 
        listed in a data frame (as produced by summarise_transcripts()).
        
        By defualt, cds sequences will be retrieved; setting cds = False
        will retrieve cDNA sequences instead."""
        
        if cds:
            seqResults = self._get_cds(list(transcriptData.TranscriptID))
        else:
            seqResults = self._get_cDNA(list(transcriptData.TranscriptID))
        
        # Convert to dictionary and check for potential issues
        resultDict = {}
        
        for result in seqResults:            
            if result['id'] != result['query']:
                raise ValueError("Sequence returned for transcript {} does not match input TranscriptID".format(result['query']))
            
            transcriptID = result['id']
            geneID = transcriptData.GeneID[transcriptData.TranscriptID == transcriptID]
            result['gene_id'] = geneID.values[0]  # Assuming only value

            resultDict[transcriptID] = result
            
        # Check if all ids returned results:
        if len(resultDict.keys()) != len(set(transcriptData.TranscriptID)):
            missingIDs = [ident for ident in transcriptData.TranscriptID if ident not in resultDict.keys()]
            
            # Check if these transcripts are simply not protein coding, while CDS seqs were requested:
            proteinCoding = self._check_protein_coding(missingIDs)
            noMatch = [ident for ident in missingIDs if proteinCoding[ident] is None]
            missingPC = [ident for ident in missingIDs if proteinCoding[ident]]
            
            if not cds and len(noMatch) == 0 and len(missingPC) == 0:
                # All missing sequences are explained
                warnings.warn('cDNA requested, but not all transcriptIDs supplied were protein coding (these were ignored)')
            else:
                if not cds and len(missingPC) != 0:
                    idtext = ', '.join(missingPC)
                    message = 'Valid protein-coding TranscriptIDs did not return cDNA results: {0}'.format(idtext)
                    warnings.warn(message)
                
                if len(noMatch) != 0:
                    idtext = ', '.join(noMatch)
                    message = 'Some TranscriptIDs could not be found: {0}'.format(idtext)
                    warnings.warn(message)
                
        return resultDict
    
    
    def get_ensemblID_from_gene_names(self, geneNames):
        accesspoint = '/lookup/symbol/homo_sapiens/'
        
        responses = self._post_request(accesspoint, geneNames, 
                                       idName = 'symbols',
                                       filters = {})
        results = []
        
        for resp in responses:
            parsed = json.loads(resp)
            
            for gene in parsed.keys():
                results.append({'GeneName': gene,
                                'GeneID': parsed[gene]['id']})
            
        return results



# --------------------------------------------------------------------------------------------------------------------------
# Other functions
# --------------------------------------------------------------------------------------------------------------------------
def seqdict_to_fasta(outname, seqdict, mode = 'wt'):
    """Write a dictionary of sequences to a fasta file: """
    with open(outname, mode = mode) as fastaFile:
        for transcriptID in seqdict.keys():
            geneID = seqdict[transcriptID]['gene_id']
            
            name = '>{0}_{1}\n'.format(geneID, transcriptID)
            seq = '{}\n'.format(seqdict[transcriptID]['seq'])
            fastaFile.writelines([name, seq])



def retry(fun, times, wait = 15, **kwargs):
    """ Attempt to run a function ('fun'), re-trying a set number of times 
    (possibly waiting a while before trying again)
    All unrecognized keyword arguments are passed to 'fun'  """
    
    if times > 1:
        try:
            result = fun(**kwargs)
        except Exception as e:
            print('{0}. Re-trying in {1}ms...'.format(type(e).__name__, wait))
            sleep(wait)
            result = retry(fun = fun, times = times - 1,
                           wait = wait, **kwargs)
            
    elif times == 1:
        # re-try one more time (this causes the correct error message to
        # be returned if this attempt fails too)
        result = fun(**kwargs)
    
    return result


# --------------------------------------------------------------------------------------------------------------------------
## Download sequences:
# --------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # --------------------------------------------------------------------------------------------------------------------------
    ## Resolve location of gene set files relative to this script:
    # --------------------------------------------------------------------------------------------------------------------------
    scriptPath = os.path.dirname(os.path.realpath(__file__))
    
    # ISG data
    ISGFile = os.path.join(scriptPath, ISG)
    ISGData = pandas.read_csv(ISGFile, index_col = 0)
    
    
    # Housekeeping gene names:
    HouseFile = os.path.join(scriptPath, HOUSEKEEPING)
    HousekeepingData = pandas.read_csv(HouseFile, sep = '\t', header = None, 
                                       names = ['GeneName', 'Accession'])
    
    
    # Change to output dir:
    outDir = os .path.join(scriptPath, DATA_PATH)
    os.chdir(outDir)
    
    
    # --------------------------------------------------------------------------------------------------------------------------
    ## Retrieve the ENSEMBL IDs for housekeeping genes:
    # --------------------------------------------------------------------------------------------------------------------------
    requester = Ensembl()
    
    housekeepingNames = list(HousekeepingData.GeneName)
    housekeepingIDs = requester.get_ensemblID_from_gene_names(housekeepingNames)
    
    housekeepingIDs = pandas.DataFrame(housekeepingIDs)
    housekeepingIDs.to_csv('HousekeepingGeneIDs.csv', index = False)
    
    
    # --------------------------------------------------------------------------------------------------------------------------
    ## Unique ENSEMBL IDs to retrieve:
    # --------------------------------------------------------------------------------------------------------------------------
    ISGnames = [name for name in ISGData.index if not name.startswith('CVR')]  # Exclude internal / newly identified genes
    
    AllNames = ISGnames + list(housekeepingIDs.GeneID)
    AllNames = list(set(AllNames))
    
    
    # --------------------------------------------------------------------------------------------------------------------------
    ## Retrieve transcripts for these genes IDs:
    # --------------------------------------------------------------------------------------------------------------------------
    
    ## Look up transcripts associated with each gene
    #   - To guard against timeouts, do this in batches:
    batchSize = 100
    writeHeader = True
    
    for start in range(0, len(AllNames), batchSize):
        stop = min(start + batchSize, len(AllNames))
        print('Fetching ids {0} to {1}'.format(start, stop))
        
        # Try up to five times to absorb brief network disconnections:
        ts = retry(requester.summarise_transcripts, times = 5, wait = 10,
                   idList = AllNames[start:stop])
        
        if writeHeader:
            ts.to_csv('TranscriptData.csv', mode = 'w', index = False, header = writeHeader)  # Create file, including header
            writeHeader = False   # Subsequent iterations should add content only
        else:
            ts.to_csv('TranscriptData.csv', mode = 'a', index = False, header = writeHeader)  # Append to existing file
    
    
    transcriptSummary = pandas.read_csv('TranscriptData.csv')
        
    
    
    ## Get coding sequences for all (canonical) transcripts:
    #   - For protein coding genes, get the cds, so we have the exact reading frame
    #   - For all other genes, get the cDNA
    #   - For genes with no canonical transcript marked, retrieve all transcripts
    #     and keep the longest one. If such a gene has both protein coding and 
    #     non protein-coding transcripts, get the longest protein-coding one.
    
    # Remove genes which have been retired from ENSEMBL - these will not have 
    #  a transcript listed:
    transcriptSummary.dropna(inplace = True)
    
    # Identify sets requiring different treatment
    canonicalProtein = pandas.DataFrame()
    canonicalOther = pandas.DataFrame()
    noncProtein = pandas.DataFrame()
    noncOther = pandas.DataFrame()
    
    for gene in set(transcriptSummary.GeneID):
        subData = transcriptSummary[transcriptSummary.GeneID == gene]
        canonData = subData[subData.Canonical == 1]
        
        if len(canonData.index) != 0:
            if len(canonData.index) > 1:
                raise ValueError('More than one canonical transcript for this gene')
            
            if all(canonData.Biotype == 'protein_coding'):
                canonicalProtein= canonicalProtein.append(canonData)
            else:
                canonicalOther = canonicalOther.append(canonData)
                
        else:
            if any(subData.Biotype == 'protein_coding'):
                protData = subData[subData.Biotype == 'protein_coding']
                noncProtein = noncProtein.append(protData)
            else:
                noncOther = noncOther.append(subData)
    
    
    # Download sequences
    open('TranscriptSequences.fasta', 'wt').close()  # Clear the file if it exists (will be appending to it below)
    
    def fetch_batch(data, cds = True, filterLength = False):
        for start in range(0, len(data), batchSize):
            stop = min(start + batchSize, len(data))
            print('Fetching ids {0} to {1}'.format(start, stop))
            
            batchData = data[start:stop]
            seqs = retry(requester.get_sequences, times = 5, wait = 10,  # Try five times to absorb brief network disconnections
                         transcriptData = batchData, cds = cds)
            
            if filterLength:
                keepSeqs = {}
                
                for gene in set(batchData.GeneID):
                    transcripts = batchData.TranscriptID[batchData.GeneID == gene]
                    longestFound = 0
                    keepID = None
                    
                    for tid in transcripts:
                        if len(seqs[tid]['seq']) > longestFound:
                            keepID = tid
                    
                    keepSeqs[keepID] = seqs[keepID]
                seqs = keepSeqs
                
            seqdict_to_fasta('TranscriptSequences.fasta', seqs, mode = 'a')
            
        return 'Done'
    
    
    fetch_batch(canonicalProtein, cds = True)
    fetch_batch(canonicalOther, cds = False)
    
    fetch_batch(noncProtein, cds = True, filterLength = True)
    fetch_batch(noncOther, cds = False, filterLength = True)
    
    
    
    # --------------------------------------------------------------------------------------------------------------------------
    ## Output:
    # --------------------------------------------------------------------------------------------------------------------------
    # Above steps already created:
    # - HousekeepingGeneIDs.csv
    # - TranscriptData.csv
    # - TranscriptSequences.csv
    #
    
    
    
    