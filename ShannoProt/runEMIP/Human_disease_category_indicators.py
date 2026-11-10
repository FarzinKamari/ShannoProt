#!/usr/bin/env python
# coding: utf-8

# # Human disease category indicators -- EMIP, etc.
# * This code integrates UniProt entries, with PICKLE PPI data, and Orphanet disease data to prepare human disease project analyses.

# In[1]:


import csv
import numpy as np
import pandas as pd
tiny = np.finfo(np.float32).tiny


# * Integrating PICKLE data with UniProt
# ##### PICKLE PPI data TXT file from:
# >* http://www.pickle.gr/Data/2.5/PICKLE2_5_UniProtNormalizedTabular-default.zip
# ##### Human proteome data TAB file from:
# >* https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&force=true&columns=id,reviewed,genes(PREFERRED),protein%20names,sequence,database(Orphanet),comment(INVOLVEMENT%20IN%20DISEASE)&compress=yes
# >###### *** Downloaded files and the code should be in the same directory

# In[2]:


pickle_ints = pd.read_csv('UniProtNormalizedTabular-default.txt', delimiter='\t', usecols=['InteractorA', 'InteractorB']) # import PICKLE data
Human_proteome = pd.read_csv('uniprot-proteome UP000005640.tab', delimiter='\t') # import human proteome data

# Converting PICKLE binary interaction format to uniprot 'interacts with' format
grp1 = pickle_ints.groupby('InteractorA')['InteractorB'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with').rename(columns={'InteractorA':'Entry'})
grp2 = pickle_ints.groupby('InteractorB')['InteractorA'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with').rename(columns={'InteractorB':'Entry'})

# Replacing interactions column of the PICKLE database
ints_in_uniprot_format = pd.concat([grp1, grp2]).groupby('Entry')['Interacts with'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with')
ints_in_uniprot_format['Interacts with'] = ints_in_uniprot_format['Interacts with'].apply(lambda x : x.split('; ')).apply(lambda x : list(set(list(x)))).apply(lambda x : '; '.join(str(elem) for elem in x))
Human_proteome = Human_proteome.merge(ints_in_uniprot_format, how='left', on='Entry')

# Generating a file for human proteome data integrated with PICKLE interactions
Human_proteome.to_csv('Human proteome with Pickle interactions.csv', index=False)


# Importing the generated CSV file and creating dictionaries and lists:
# * indx is a dictionary that will be used to find a given protein entry and return its index.
# * ints is the list of interaction lists (list of lists).
# * seqs is the list for sequence strings.

# In[3]:


indx = {} # entry indexes dictionary
ints = [] # list of interactions
n_ints =[] # list of number of interactions
seqs = [] # list of sequence strings
Entry= [] # list of Entries
prot_name = [] # list of protein names
orpha_uniprot = [] # list of orphanet numbers
involvement = [] # list of involvement in diseases
status = [] # list of Sprot/trembl status
gene_name = [] # list of Gene names (primary)

with open('Human proteome with Pickle interactions.csv', newline = '') as csvfile:
    n = 0 # number of row
    for row in csv.DictReader(csvfile):
        indx[row['Entry']] = n
        # create list of lists for interactions and appending other lists as needed:
        ints.append([entry for entry in row['Interacts with'].split('; ')] if row['Interacts with'] != '' else [])
        n_ints.append(row['Interacts with'].count('; ') + 1 if row['Interacts with'] != '' else 0)
        seqs.append(row['Sequence'])
        Entry.append(row['Entry'])
        prot_name.append(row['Protein names'])
        orpha_uniprot.append(row['Cross-reference (Orphanet)'])
        involvement.append(row['Involvement in disease'])
        status.append(row['Status'])
        gene_name.append(row['Gene names  (primary )'])
        n += 1


# All alphabet letters are counted in sequence strings forming the rep_ent matrix whose rows are representative of protein entries and columns are representative of A-Z letters.

# In[4]:


rep_ent = np.zeros((n, 22), dtype = np.float32) # matrix of number of repitition of amino acids in entry sequences initialized with zeros
for i in range(n):
    rep_ent[(i,  0)] = seqs[i].count('A') # count A string in Sequence
    rep_ent[(i,  1)] = seqs[i].count('C') # count C string in Sequence
    rep_ent[(i,  2)] = seqs[i].count('D') + seqs[i].count('B') / 2 # count D string in Sequence
    #(Since B residue is either D or N, half of its frequency has been added here)
    rep_ent[(i,  3)] = seqs[i].count('E') + seqs[i].count('Z') / 2 # count E string in Sequence
    #(Since Z residue is either E or Q, half of its frequency has been added here)
    rep_ent[(i,  4)] = seqs[i].count('F') # count F string in Sequence
    rep_ent[(i,  5)] = seqs[i].count('G') # count G string in Sequence
    rep_ent[(i,  6)] = seqs[i].count('H') # count H string in Sequence
    rep_ent[(i,  7)] = seqs[i].count('I') + seqs[i].count('J') / 2 # count I string in Sequence
    #(Since J residue is either I or L, half of its frequency has been added here)
    rep_ent[(i, 8)] = seqs[i].count('K') # count K string in Sequence
    rep_ent[(i, 9)] = seqs[i].count('L') + seqs[i].count('J') / 2 # count L string in Sequence
    #(Since J residue is either I or L, half of its frequency has been added here)
    rep_ent[(i, 10)] = seqs[i].count('M') # count M string in Sequence
    rep_ent[(i, 11)] = seqs[i].count('N') + seqs[i].count('B') / 2 # count N string in Sequence
    #(Since B residue is either D or N, half of its frequency has been added here)
    rep_ent[(i, 12)] = seqs[i].count('O') # count O string in Sequence
    rep_ent[(i, 13)] = seqs[i].count('P') # count P string in Sequence
    rep_ent[(i, 14)] = seqs[i].count('Q') + seqs[i].count('Z') / 2 # count Q string in Sequence
    #(Since Z residue is either E or Q, half of its frequency has been added here)
    rep_ent[(i, 15)] = seqs[i].count('R') # count R string in Sequence
    rep_ent[(i, 16)] = seqs[i].count('S') # count S string in Sequence
    rep_ent[(i, 17)] = seqs[i].count('T') # count T string in Sequence
    rep_ent[(i, 18)] = seqs[i].count('U') # count U string in Sequence
    rep_ent[(i, 19)] = seqs[i].count('V') # count V string in Sequence
    rep_ent[(i, 20)] = seqs[i].count('W') # count W string in Sequence
    rep_ent[(i, 21)] = seqs[i].count('Y') # count Y string in Sequence


# Calculating disease indicators: CAIR and EMIP

# In[5]:


len_ent = rep_ent @ np.ones((22, 1), dtype = np.float32) # the length of Entry sequences is calculated by ones-matrix multiplication
P = rep_ent / (len_ent) # probability of each amino acid for Entries
cair_ent = -np.sum(P * (np.log2(P + tiny) / np.log2(22)), axis = 1, keepdims = True) # calculating CAIRs of Entries (tiny is used to avoid log0)
pi_int = np.zeros((n, 1), dtype = np.float32) # creating zeros matrix
rep_net = np.copy(rep_ent) # hard copying rep_ent for rep_net (i.e. frequency of residues in a PPI network) calculation
for i in range(n):
    for entry in ints[i]:
        try:
            ind = indx[entry]
            pi_int[(i, 0)] += cair_ent[(ind, 0)] * len_ent[(ind, 0)] # pi_int matrix is filled in by summation of interactors' PIs
            rep_net[(i), 0:22] += rep_ent[(ind), 0:22] # rep_net is filled in by adding residue frequencies of interactors 
        except KeyError:
            None
len_net = rep_net @ np.ones((22, 1), dtype = np.float32) # calculating accumulative residue length of networks (by ones-matrix multiplication)
P = rep_net / (len_net) # probability of each amino acid for networks
cair_net = -np.sum(P * (np.log2(P + tiny) / np.log2(22)), axis = 1, keepdims = True) # calculating CAIRs of networks (tiny is used to avoid log0)
net_inf = len_net * cair_net # network information
emip = net_inf - pi_int # EMIP


# * Integrating the data with Orphanet disease database and gene age data which are available open-source at:
# ##### Orphanet database XML file from:
# >* http://www.orphadata.org/data/xml/en_product9_prev.xml
# >###### *** This file should be saved as a CSV format file using Excel.
# ##### Gene age consensus article:
# >* Liebeskind BJ, McWhite CD, Marcotte EM. Towards consensus gene ages. Genome biology and evolution. 2016 Jun 1;8(6):1812-23.
# >* https://github.com/marcottelab/Gene-Ages/blob/master/Main/main_HUMAN.csv
# 
# >###### *** Downloaded files and the code should be in the same directory

# In[6]:


Age = pd.read_csv('main_HUMAN.csv', usecols=['Entry', 'modeAge']) # adding gene ages

## preparing and adding the Orphanet data as described in Methods section
# preparing Orphanet raw (unfiltered) data
orphadata_unprocessed = pd.read_csv('en_product9_prev.csv', encoding = "ISO-8859-1",
                                    usecols=['OrphaNumber', 'ExpertLink', 'Name', 'Name4', 'Source', 'Name12', 'Name18', 'Name21',
                                             'Name24']).rename(columns={'OrphaNumber': 'OrphaNumber (source: Orphanet)',
                                                                        'Name': 'Disease name (source: Orphanet)', 'Name4': 'Disease type',
                                                                        'Name12': 'Occurence type', 'Name18': 'Occurence value',
                                                                        'Name21': 'Geo distrib', 'Name24': 'Validity'}) # reading columns renaming them
orphadata_unprocessed['OrphaNumber (source: Orphanet)'] = orphadata_unprocessed['OrphaNumber (source: Orphanet)'].astype(str) # raw data of Orphanet

# preparing Orphanet processed (filtered) data
orphadata = pd.read_csv('en_product9_prev.csv', encoding = "ISO-8859-1",
                        usecols=['OrphaNumber', 'ExpertLink', 'Name', 'Name4', 'Source', 'Name12', 'Name18', 'Name21',
                                 'Name24']).rename(columns={'OrphaNumber': 'OrphaNumber (source: Orphanet)',
                                                            'Name': 'Disease name (source: Orphanet)', 'Name4': 'Disease type',
                                                            'Name12': 'Occurence type', 'Name18': 'Occurence value',
                                                            'Name21': 'Geo distrib', 'Name24': 'Validity'}) # reading columns again renaming them 
# applying filters as explained in Methods section
orphadata = orphadata[orphadata['Geo distrib'] == 'Worldwide'] # filtering out all locally distributed diseases -- Only keeps 'Worldwide'
orphadata = orphadata[orphadata['Occurence value'] != 'Unknown'] # filtering out 'Unknown' occurences
orphadata = orphadata[orphadata['Occurence value'] != 'Not yet documented'] # filtering out missing occurence values
orphadata = orphadata[orphadata['Occurence type'] != 'Lifetime Prevalence'] # filtering out 'Lifetime Prevalence'
orphadata = orphadata[orphadata['Occurence type'] != 'Cases/families'] # filtering out 'Cases/families'
orphadata = orphadata.dropna(subset=['Occurence value']) # removing NAs in Occurence column
orphadata = orphadata.drop_duplicates('OrphaNumber (source: Orphanet)', keep='first') # preference is implied as below: (see methods section)
                                                                                        # 1- Incidence, 2- Prevalence at birth, 3- Point prevalence
orphadata['OrphaNumber (source: Orphanet)'] = orphadata['OrphaNumber (source: Orphanet)'].astype(str)


# ### Getting output files

# In[7]:


## the first output (out1) file only designates the 'disease probability' according to the UniProt 'Involvement in disease' column
prim_out = pd.DataFrame({'Status': status,
                        'Entry': Entry,
                        'Gene name': gene_name,
                        'Protein name': prot_name,
                        'Protein length': len_ent[:, 0],
                        'CAIR (in pits)': cair_ent[:, 0],
                        'Number of interactions': n_ints,
                        'EMIP (in pits)': emip[:, 0],
                        'Orpha no. (source: UniProt)': orpha_uniprot,
                        'Involvement in disease (source: UniProt)': involvement}) # creating the primary output dataframe
prim_out = prim_out.merge(Age, how='left', on='Entry').rename(columns={'modeAge': 'Gene age'}) # adding gene age column
out1 = prim_out.copy() # hard-copying to write out1
out1.to_csv('w%entries_w%diseases_w%uniprot_wo%orpha.csv', index = False) # writing the CSV out1

## the second output (out2) file merges UNPROCESSED Orphadata separating the data with respect to Orpha numbers; represents all of the available disease data
prim_out['Orpha no. (source: UniProt)'] = prim_out['Orpha no. (source: UniProt)'].apply(lambda x : x.split(';')[:-1]).apply(list) # preparing for the 'explode' function
prim_out = prim_out.explode('Orpha no. (source: UniProt)') # requires pandas version 0.25 or later; separates the data with respect to Orpha numbers
out2 = prim_out.merge(orphadata_unprocessed, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)')
out2 = out2.reindex(columns=['Orpha no. (source: UniProt)',
                             'Involvement in disease (source: UniProt)', 
                             'OrphaNumber (source: Orphanet)',
                             'Disease name (source: Orphanet)',
                             'Disease type', 
                             'Occurence value',
                             'Occurence type',
                             'Geo distrib',
                             'Validity',
                             'Source',
                             'Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'CAIR (in pits)',
                             'Number of interactions',
                             'EMIP (in pits)',
                             'Gene age',
                             'ExpertLink']) # reindexing columns
out2['Orpha no. (source: UniProt)'] = out2['Orpha no. (source: UniProt)'].astype(float) # retyping Orpha numbers
out2 = out2.sort_values(by=['Orpha no. (source: UniProt)']) # sorting based on Orpha numbers
out2.to_csv('w%entries_w%diseases(expanded)_w%uniprot_w%orpha(unprocessed).csv', index = False) # writing the CSV output

## the third output (out3) file merges PROCESSED Orphadata deleting all proteins w/o an Orpha number
out3 = prim_out.merge(orphadata, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)') # merges Orpha numbers from two sources
out3 = out3.dropna(subset=['OrphaNumber (source: Orphanet)']) # dropping NAs
out3 = out3.reindex(columns=['OrphaNumber (source: Orphanet)',
                             'Disease name (source: Orphanet)',
                             'Involvement in disease (source: UniProt)',
                             'Disease type',
                             'Occurence value',
                             'Occurence type',
                             'Geo distrib',
                             'Validity',
                             'Source',
                             'Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'CAIR (in pits)',
                             'Number of interactions',
                             'EMIP (in pits)',
                             'Gene age',
                             'ExpertLink']) # reindexing columns
out3.to_csv('w%entries(dropped)_w%diseases(expanded)_w%uniprot_w%orpha(processed).csv', index=False) # writing the CSV output

## the forth output (out4) file merges PROCESSED Orphadata and calculates total occurence for all proteins; out4 is the final output
out4 = prim_out.merge(orphadata, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)') # merges Orpha numbers
out4['Orpha no. (source: UniProt)'] = out4['Orpha no. (source: UniProt)'].fillna('0') # filling NAs with zeros to avoid unequality of NAs
out4['OrphaNumber (source: Orphanet)'] = out4['OrphaNumber (source: Orphanet)'].fillna('0') # filling NAs with zeros to avoid unequality of NAs
out4['Occurence value'] = out4.apply(lambda row : row['Occurence value'] if row['Orpha no. (source: UniProt)'] == row['OrphaNumber (source: Orphanet)'] else float(tiny), axis=1) # separating out NAs (tiny as a Python float, so that the occurrences are summed in double precision)
out4['Occurence value'] = out4['Occurence value'].fillna(0) # Giving zeros for the practicability of summation
out4['Occurence value'] = out4['Occurence value'].replace({'<1 / 1 000 000' : 0.000001, '1-9 / 100 000': 0.00005, '1-9 / 1 000 000': 0.000005,
                                                           '1-5 / 10 000': 0.0003, '>1 / 1000': 0.001, '6-9 / 10 000': 0.00075}) # quantifying the ranked Orpha data with means of the rank

Tot_occur = out4[out4['Occurence value'] != 0].groupby('Entry')['Occurence value'].apply(sum).reset_index(name='Total occurence value') # summing up total occurences
Tot_occur['Total occurence value'] = Tot_occur['Total occurence value'].apply(lambda x : 'NA' if x < 0.000000001 else x) # returning NA values back
out4 = out1.merge(Tot_occur, how='left', on='Entry') # merging total occurences
out4['Total occurence value'] = out4['Total occurence value'].fillna(0) # Giving zeros to non-diseases to distinguish it from real NAs
out4 = out4.reindex(columns=['Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'CAIR (in pits)',
                             'Number of interactions',
                             'EMIP (in pits)',
                             'Orpha no. (source: UniProt)',
                             'Total occurence value',
                             'Involvement in disease (source: UniProt)',
                             'Gene age']) # reindexing columns
out4.to_csv('w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv', index=False) # writing the CSV output


# In[ ]:





# In[ ]:




