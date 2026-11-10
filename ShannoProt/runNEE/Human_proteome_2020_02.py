#!/usr/bin/env python
# coding: utf-8

# # Human proteome of UniProt release 2020_02 -- input table for NEE (formerly EMIP)
# * In 2020 the human proteome table (UP000005640) was downloaded from the UniProt website (see runEMIP). The website now serves only the current release, so this module rebuilds the reviewed (Swiss-Prot) part of that table from the archived release 2020_02.
# * PICKLE 2.5 interactions are among reviewed entries, so the Swiss-Prot part is all that the NEE calculations need.
# * The output TAB file has the same columns as the 2020 download (plus the secondary accessions) and can be read by Human_disease_category_indicators.py.
# ##### UniProt release 2020_02 (Swiss-Prot) from:
# >* https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2020_02/knowledgebase/uniprot_sprot-only2020_02.tar.gz
# >###### *** uniprot_sprot.dat.gz should be extracted from the archive (see README)

# In[1]:


import gzip
import re
import pandas as pd


# In[2]:


def dat_to_proteome_table(infile='uniprot_sprot.dat.gz', outfile='Human proteome 2020_02.tab', taxid='9606'): # refer to runNEE
    rows = [] # one dictionary per human entry
    entry = {} # the entry being read
    with gzip.open(infile, 'rt', encoding='utf-8') as fileread:
        for line in fileread:
            code = line[:2] # line code of the UniProt DAT format
            if code == 'ID':
                entry = {'acc': [], 'gene': '', 'name': '', 'seq': [], 'orpha': [], 'disease': [], 'proteomes': [], 'geneid': [], 'ensembl': [], 'ox': '', 'status': line.split()[2].rstrip(';')}
                in_disease = False # True while reading a "DISEASE" comment
                in_seq = False # True while reading the sequence block
            elif code == 'AC':
                entry['acc'] += [a.strip() for a in line[5:].split(';') if a.strip()] # the first accession is the primary one
            elif code == 'OX':
                entry['ox'] = re.search(r'NCBI_TaxID=(\d+)', line).group(1) # organism ID
            elif code == 'DE' and line[5:].startswith('RecName: Full=') and not entry['name']:
                entry['name'] = re.sub(r' \{ECO:.*?\}', '', line[19:].strip()).rstrip(';') # recommended protein name w/o evidence tags
            elif code == 'GN' and 'Name=' in line and not entry['gene']:
                entry['gene'] = re.sub(r' \{ECO:.*?\}', '', line.split('Name=')[1].split(';')[0]).strip() # primary gene name
            elif code == 'DR' and line[5:].startswith('Orphanet;'):
                entry['orpha'].append(line[5:].split(';')[1].strip()) # Orpha numbers cross-referenced by UniProt
            elif code == 'DR' and line[5:].startswith('GeneID;'):
                entry['geneid'].append(line[5:].split(';')[1].strip()) # NCBI GeneID (used for the PubMed counts)
            elif code == 'DR' and line[5:].startswith('Ensembl;'):
                ens = line[5:].split(';')[3].strip().split('.')[0] # Ensembl gene ID (used for HuRI)
                if ens not in entry['ensembl']:
                    entry['ensembl'].append(ens)
            elif code == 'DR' and line[5:].startswith('Proteomes;'):
                entry['proteomes'].append(line[5:].split(';')[1].strip()) # e.g. UP000005640 (the human reference proteome)
            elif code == 'CC':
                if line[5:].startswith('-!- DISEASE:'):
                    in_disease = True
                    entry['disease'].append(line[9:].strip())
                elif line[5:].startswith('-!-') or line[5:].startswith('---'):
                    in_disease = False
                elif in_disease:
                    entry['disease'][-1] += ' ' + line[5:].strip() # continuation lines of the DISEASE comment
            elif code == 'SQ':
                in_seq = True
            elif code == '//':
                if entry['ox'] == taxid:
                    rows.append({'Entry': entry['acc'][0],
                                 'Status': entry['status'].lower(),
                                 'Gene names  (primary )': entry['gene'],
                                 'Protein names': entry['name'],
                                 'Sequence': ''.join(entry['seq']),
                                 'Cross-reference (Orphanet)': ''.join(o + ';' for o in entry['orpha']), # same format as the 2020 download (e.g. "93;1234;")
                                 'Involvement in disease': ' '.join(entry['disease']),
                                 'Secondary accessions': '; '.join(entry['acc'][1:]), # used to reconcile obsolete identifiers
                                 'Proteomes': '; '.join(entry['proteomes']), # the 2020 table contained the entries of UP000005640 only
                                 'GeneID': '; '.join(entry['geneid']),
                                 'Ensembl genes': '; '.join(entry['ensembl'])})
                in_seq = False
            elif in_seq:
                entry['seq'].append(line.replace(' ', '').strip()) # sequence lines
    out = pd.DataFrame(rows) # creating the output dataframe
    out.to_csv(outfile, sep='\t', index=False) # writing the TAB file
    print(len(out), 'entries of organism', taxid, 'written to', outfile)
    return out

