#!/usr/bin/env python
# coding: utf-8

# # NEE (formerly EMIP) and its decomposition -- human proteome
# * This module recalculates NEE as the released runEMIP code does (float32), checks the values against Table S3 of the preprint, and adds:
# >* L_N and CEI_N: length and CEI of a protein pooled with its interaction partners
# >* the own term L*CEI and the network term M = NEE - L*CEI = L_N*D_JS, where D_JS is the length-weighted Jensen-Shannon divergence of the compositions in the neighborhood
# >* CEI_MM and NEE_MM: Miller-Madow bias correction (as described in Methods section)
# * NEE = L_N*CEI_N - (sum of L_j*CEI_j over the partners j), in base-22 units
# * Interactions are read as in runEMIP: all PICKLE partners count for the number of interactions (self-interactions included); only partners present in the proteome table are pooled.
# ##### Input files:
# >* Human proteome 2020_02.tab (from Human_proteome_2020_02.py)
# >* UniProtNormalizedTabular-default.txt (PICKLE 2.5) from http://www.pickle.gr/Data/2.5/PICKLE2_5_UniProtNormalizedTabular-default.zip
# >* Table S3 of the preprint (w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv in Supplementary_materials; only for the check)
# ##### Output file:
# >* NEE decomposition.csv

# In[1]:


import numpy as np
import pandas as pd
from scipy import sparse
tiny = np.finfo(np.float32).tiny
residues = 'ACDEFGHIKLMNOPQRSTUVWY' # column order of the residue matrix (as in runEMIP)


# In[2]:


def residue_counts(seqs, dtype=np.float64): # refer to nee_decomposition
    rep = np.zeros((len(seqs), 22), dtype = dtype) # matrix of number of repetition of residues (rows: entries, columns: residues)
    for i, seq in enumerate(seqs):
        for j, aa in enumerate(residues):
            rep[i, j] = seq.count(aa) # count each residue in the sequence
        rep[i, 2] += seq.count('B') / 2 # D (since B residue is either D or N, half of its frequency is added to each)
        rep[i, 11] += seq.count('B') / 2 # N
        rep[i, 3] += seq.count('Z') / 2 # E (since Z residue is either E or Q, half of its frequency is added to each)
        rep[i, 14] += seq.count('Z') / 2 # Q
        rep[i, 7] += seq.count('J') / 2 # I (since J residue is either I or L, half of its frequency is added to each)
        rep[i, 9] += seq.count('J') / 2 # L
    return rep


def cei_of(rep): # plug-in CEI of each row of a residue matrix (float64, no tiny needed)
    length = rep.sum(axis = 1, keepdims = True) # lengths
    P = np.divide(rep, length, out = np.zeros_like(rep), where = length > 0) # residue probabilities
    logP = np.log(P, out = np.zeros_like(P), where = P > 0) # 0*log0 is taken as 0
    return -np.sum(P * logP, axis = 1) / np.log(22) # Shannon entropy divided by ln22, i.e. CEI


def read_interactions(pickle_file): # partners of each entry as in runEMIP ('Interacts with' format)
    pickle_ints = pd.read_csv(pickle_file, delimiter='\t', usecols=['InteractorA', 'InteractorB']) # import PICKLE data
    ints = {} # entry -> set of partners
    for a, b in zip(pickle_ints.InteractorA, pickle_ints.InteractorB):
        ints.setdefault(a, set()).add(b) # both directions, as the grouping of runEMIP does
        ints.setdefault(b, set()).add(a) # (a self-interaction adds the entry to its own partners)
    return ints


def released_emip(rep, ints_list, indx): # the released algorithm of runEMIP, in float32 (used for the check against Table S3 of the preprint)
    rep_ent = rep.astype(np.float32) # matrix of residue counts
    len_ent = rep_ent @ np.ones((22, 1), dtype = np.float32) # lengths by ones-matrix multiplication
    P = rep_ent / (len_ent) # probability of each amino acid for entries
    cair_ent = -np.sum(P * (np.log2(P + tiny) / np.log2(22)), axis = 1, keepdims = True) # CAIRs (CEIs) of entries
    pi_int = np.zeros((len(rep_ent), 1), dtype = np.float32) # summation of partners' L*CEI
    rep_net = np.copy(rep_ent) # residue counts of the networks
    for i in range(len(rep_ent)):
        for entry in ints_list[i]:
            try:
                ind = indx[entry]
                pi_int[(i, 0)] += cair_ent[(ind, 0)] * len_ent[(ind, 0)]
                rep_net[(i), 0:22] += rep_ent[(ind), 0:22]
            except KeyError:
                None
    len_net = rep_net @ np.ones((22, 1), dtype = np.float32) # network lengths
    P = rep_net / (len_net) # probability of each amino acid for networks
    cair_net = -np.sum(P * (np.log2(P + tiny) / np.log2(22)), axis = 1, keepdims = True) # CAIRs of networks
    return (len_net * cair_net - pi_int)[:, 0], cair_ent[:, 0] # EMIP and CAIR as in runEMIP


# In[3]:


def nee_decomposition(proteome_file='Human proteome 2020_02.tab', pickle_file='UniProtNormalizedTabular-default.txt',
                      outfile='NEE decomposition.csv', proteome_id='UP000005640'): # refer to runNEE
    prot = pd.read_csv(proteome_file, sep='\t', keep_default_na=False) # human Swiss-Prot entries of release 2020_02
    prot = prot[prot['Proteomes'].str.contains(proteome_id)].reset_index(drop=True) # only the entries of the reference proteome (as in 2020)
    indx = {e: i for i, e in enumerate(prot.Entry)} # entry indexes dictionary
    ints = read_interactions(pickle_file) # PICKLE partners
    ints_list = [sorted(ints.get(e, set())) for e in prot.Entry] # list of interactions for each entry
    rep = residue_counts(prot.Sequence.tolist()) # residue counts (float64)

    # adjacency matrix of the partners that are present in the proteome table (these are pooled)
    rows = [i for i, part in enumerate(ints_list) for p in part if p in indx]
    cols = [indx[p] for part in ints_list for p in part if p in indx]
    adj = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape = (len(prot), len(prot)))

    length = rep.sum(axis = 1) # protein length L
    cei = cei_of(rep) # CEI of each protein
    m_obs = (rep > 0).sum(axis = 1) # number of residue types observed in each protein
    rep_net = rep + adj @ rep # residues of the protein pooled with its partners
    len_net = rep_net.sum(axis = 1) # L_N
    cei_net = cei_of(rep_net) # CEI_N
    m_net = (rep_net > 0).sum(axis = 1) # residue types observed in the pooled neighborhood
    own = length * cei # own term L*CEI (= S/log2(22), total compositional entropy in base-22 units)
    partners = adj @ own # summation of the partners' L*CEI
    nee = len_net * cei_net - partners # NEE
    emip_f32, cair_f32 = released_emip(rep, ints_list, indx) # released algorithm, for the check

    out = pd.DataFrame({'Entry': prot.Entry,
                        'Protein length': length,
                        'CEI': cei,
                        'CEI_MM': cei + (m_obs - 1) / (2 * length * np.log(22)), # Miller-Madow corrected CEI
                        'Residue types observed': m_obs,
                        'Number of interactions': [len(part) for part in ints_list], # all distinct PICKLE partners (as in runEMIP)
                        'Partners in proteome': np.asarray(adj.sum(axis = 1))[:, 0], # partners that are pooled
                        'Self-interaction': [e in ints.get(e, set()) for e in prot.Entry],
                        'L_N': len_net,
                        'CEI_N': cei_net,
                        'Own term (L*CEI)': own,
                        'Network term M': nee - own, # M = L_N*D_JS >= 0
                        'D_JS': (nee - own) / len_net, # length-weighted Jensen-Shannon divergence
                        'NEE': nee,
                        'NEE_MM': nee + (m_net - 1) / (2 * np.log(22)) - adj @ ((m_obs - 1) / (2 * np.log(22))), # Miller-Madow corrected NEE
                        'NEE (released float32 code)': emip_f32,
                        'CEI (released float32 code)': cair_f32}) # creating the output dataframe
    out.to_csv(outfile, index=False) # writing the CSV output
    print(len(out), 'entries;', int(adj.sum()), 'pooled partner links written to', outfile)
    return out


def check_against_table_s3(nee_file='NEE decomposition.csv', table_s3='../Supplementary_materials/EMIP_supplementary_files/w%entries_w%diseases(accumulated)_w%uniprot_w%orpha(processed).csv'): # reproduces the values of 2020 before anything else is trusted
    from scipy import stats
    out = pd.read_csv(nee_file) # recalculated values
    s3 = pd.read_excel(table_s3) if table_s3.endswith('.xlsx') else pd.read_csv(table_s3, low_memory=False) # Table S3 of the preprint (CSV, or the same table as XLSX)
    s3 = s3[s3['Status'] == 'reviewed'][['Entry', 'Protein length', 'CAIR (in pits)', 'Number of interactions', 'EMIP (in pits)']]
    chk = out.merge(s3, on='Entry', suffixes=('', ' (S3)')) # entries of both tables
    rel = lambda a, b: np.abs(a - b) / np.abs(b).clip(lower = 1e-12) # relative difference
    res = {'entries compared': len(chk),
           'S3 reviewed entries': len(s3),
           'length agreement': float((chk['Protein length'] == chk['Protein length (S3)']).mean()),
           'degree agreement': float((chk['Number of interactions'] == chk['Number of interactions (S3)']).mean()),
           'Spearman NEE (float32) vs S3 EMIP': float(stats.spearmanr(chk['NEE (released float32 code)'], chk['EMIP (in pits)'])[0]),
           'max rel. diff NEE (float32) vs S3 EMIP': float(rel(chk['NEE (released float32 code)'], chk['EMIP (in pits)']).max()),
           'Spearman NEE (float64) vs S3 EMIP': float(stats.spearmanr(chk['NEE'], chk['EMIP (in pits)'])[0]),
           'median rel. diff NEE (float64) vs S3 EMIP': float(rel(chk['NEE'], chk['EMIP (in pits)']).median()),
           'max rel. diff NEE (float64) vs S3 EMIP': float(rel(chk['NEE'], chk['EMIP (in pits)']).max()),
           'max abs. diff CEI (float64) vs S3 CAIR': float(np.abs(chk['CEI'] - chk['CAIR (in pits)']).max())}
    for k, v in res.items():
        print(k + ':', v)
    return res, chk

