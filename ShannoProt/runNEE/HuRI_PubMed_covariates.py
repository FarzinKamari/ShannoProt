#!/usr/bin/env python
# coding: utf-8

# # Study-bias covariates of human proteins -- HuRI degree, NEE on HuRI and PubMed counts
# * HuRI (Luck et al. 2020) is a systematic binary interactome map and therefore free of the literature bias of curated interactomes such as PICKLE. NEE is recalculated on the HuRI network with the same formula as in NEE_decomposition.py (self-interactions included, as in runEMIP).
# * PubMed counts: number of distinct PubMed articles linked to the gene of each protein (NCBI gene2pubmed), (a) for articles created up to 31 December 2019 (PICKLE 2.5 was released in December 2019) and (b) for all articles at the download date.
# * Focused PubMed counts: the same up to 2019, but only articles linked to at most 10 (or 100) human genes. Large interaction screens are linked to thousands of genes each (e.g. Rolland et al. 2014, PMID 25416956: 4,399 human genes), so that the plain count also grows with the number of reported interactions.
# * Ensembl gene IDs (HuRI) and NCBI GeneIDs (gene2pubmed) are mapped to UniProt accessions with the cross-references of UniProt release 2020_02 itself (DR lines, see Human_proteome_2020_02.py), i.e. the release of the analysis.
# ##### Input files:
# >* NEE decomposition.csv and Human proteome 2020_02.tab (from NEE_decomposition.py and Human_proteome_2020_02.py)
# >* HuRI.tsv from http://www.interactome-atlas.org/download
# >* gene2pubmed.gz from https://ftp.ncbi.nlm.nih.gov/gene/DATA/
# ##### Output file:
# >* Study bias covariates.csv

# In[1]:


import numpy as np
import pandas as pd
from scipy import sparse
from NEE_decomposition import residue_counts, cei_of
pmid_2019 = 31887851 # the last PubMed record created on 31 December 2019 (E-utilities, [crdt])
max_genes = [10, 100] # focused articles: linked to at most this many human genes


# In[2]:


def cross_references(prot, column): # UniProt accession -> external IDs, from the DR lines of the release
    pairs = [(e, i.strip()) for e, ids in zip(prot.Entry, prot[column]) for i in ids.split(';') if i.strip()]
    return pd.DataFrame(pairs, columns=['Entry', 'ID']).drop_duplicates()


def huri_partners(huri_file, ens_map): # partners of each entry in HuRI (both directions)
    huri = pd.read_csv(huri_file, sep='\t', header=None, names=['A', 'B']) # Ensembl gene pairs (no header)
    huri = huri.merge(ens_map.rename(columns={'Entry': 'Entry_A', 'ID': 'A'}), on='A') # mapping both sides to UniProt
    huri = huri.merge(ens_map.rename(columns={'Entry': 'Entry_B', 'ID': 'B'}), on='B')
    partners = {}
    for a, b in zip(huri.Entry_A, huri.Entry_B):
        partners.setdefault(a, set()).add(b) # a homodimer adds the entry to its own partners
        partners.setdefault(b, set()).add(a)
    return partners, huri


# In[3]:


def study_bias_covariates(nee_file='NEE decomposition.csv', proteome_file='Human proteome 2020_02.tab', huri_file='HuRI.tsv',
                          gene2pubmed_file='gene2pubmed.gz',
                          outfile='Study bias covariates.csv'): # refer to runNEE
    nee = pd.read_csv(nee_file) # entries of the analysis
    prot = pd.read_csv(proteome_file, sep='\t', keep_default_na=False).set_index('Entry').loc[nee.Entry].reset_index() # same order
    indx = {e: i for i, e in enumerate(nee.Entry)} # entry indexes dictionary

    ## HuRI degree and NEE on HuRI
    ens_map = cross_references(prot, 'Ensembl genes') # Entry -> Ensembl gene
    partners, huri = huri_partners(huri_file, ens_map)
    rows = [indx[a] for a in partners for b in partners[a]]
    cols = [indx[b] for a in partners for b in partners[a]]
    adj = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape = (len(nee), len(nee))) # HuRI adjacency matrix
    rep = residue_counts(prot.Sequence.tolist()) # residue counts
    own = rep.sum(axis = 1) * cei_of(rep) # L*CEI of each protein
    rep_net = rep + adj @ rep # residues of the protein pooled with its HuRI partners
    nee_huri = rep_net.sum(axis = 1) * cei_of(rep_net) - adj @ own # NEE on HuRI

    ## PubMed counts
    gene_map = cross_references(prot, 'GeneID') # Entry -> NCBI GeneID
    g2p = pd.read_csv(gene2pubmed_file, sep='\t', dtype=str) # columns: #tax_id, GeneID, PubMed_ID
    g2p = g2p[g2p['#tax_id'] == '9606'].drop(columns='#tax_id').drop_duplicates() # human genes only
    genes_per_article = g2p.groupby('PubMed_ID').GeneID.nunique() # human genes linked to each article
    g2p = g2p.merge(gene_map.rename(columns={'ID': 'GeneID'}), on='GeneID') # adding the UniProt entries
    g2p['Genes in article'] = g2p.PubMed_ID.map(genes_per_article)
    g2p['PubMed_ID'] = g2p.PubMed_ID.astype(int)
    pubmed_all = g2p.groupby('Entry').PubMed_ID.nunique() # distinct articles per entry (several genes of an entry are pooled)
    pubmed_2019 = g2p[g2p.PubMed_ID <= pmid_2019].groupby('Entry').PubMed_ID.nunique()
    focused = {k: g2p[(g2p.PubMed_ID <= pmid_2019) & (g2p['Genes in article'] <= k)].groupby('Entry').PubMed_ID.nunique()
               for k in max_genes} # focused articles up to 2019

    out = pd.DataFrame({'Entry': nee.Entry,
                        'Ensembl gene mapped': nee.Entry.isin(set(ens_map.Entry)),
                        'HuRI degree': [len(partners.get(e, ())) for e in nee.Entry],
                        'NEE on HuRI': nee_huri,
                        'GeneID mapped': nee.Entry.isin(set(gene_map.Entry)),
                        'PubMed articles to 2019': nee.Entry.map(pubmed_2019).fillna(0).astype(int),
                        'PubMed articles (all)': nee.Entry.map(pubmed_all).fillna(0).astype(int),
                        'PubMed focused articles to 2019': nee.Entry.map(focused[10]).fillna(0).astype(int),
                        'PubMed articles to 2019 on at most 100 genes': nee.Entry.map(focused[100]).fillna(0).astype(int)}) # creating the output dataframe
    out.to_csv(outfile, index=False) # writing the CSV output
    print('HuRI pairs mapped:', len(huri), '; entries with a HuRI partner:', int((out['HuRI degree'] > 0).sum()),
          '; entries with PubMed articles to 2019:', int((out['PubMed articles to 2019'] > 0).sum()))
    return out

