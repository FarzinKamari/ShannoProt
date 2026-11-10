#!/usr/bin/env python
# coding: utf-8

# # Null models of NEE (formerly EMIP) -- composition shuffle and degree-preserving rewiring
# * Composition shuffle: each protein receives the residue frequencies of another protein (a random permutation of all proteins), scaled to its own length; lengths and interactions are kept. If NEE separated the disease groups through amino-acid composition, the shuffle would weaken the separation.
# * Degree-preserving rewiring: the interactions between proteins are rewired by random edge switches (igraph "rewire"; every protein keeps its number of partners, self-interactions are kept); compositions are kept. If the identity of the partners mattered, rewiring would weaken the separation.
# * For each replicate, NEE is recalculated and its AUC (any rare disease vs none) and Kendall's tau_b (ordered categories) are recorded.
# ##### Input files:
# >* Human proteome 2020_02.tab, UniProtNormalizedTabular-default.txt (PICKLE 2.5), Human analysis data.csv (from Human_indicator_statistics.prepare_data)
# ##### Output file:
# >* NEE null models.csv

# In[1]:


import random
import numpy as np
import pandas as pd
import igraph
from scipy import sparse, stats
from sklearn.metrics import roc_auc_score
from NEE_decomposition import residue_counts, cei_of, read_interactions


# In[2]:


def nee_of(rep, adj): # NEE of every protein for residue counts rep and partner matrix adj
    own = rep.sum(axis = 1) * cei_of(rep) # L*CEI
    rep_net = rep + adj @ rep # pooled neighborhoods
    return rep_net.sum(axis = 1) * cei_of(rep_net) - adj @ own


def nee_null_models(proteome_file='Human proteome 2020_02.tab', pickle_file='UniProtNormalizedTabular-default.txt',
                    data_file='Human analysis data.csv', outfile='NEE null models.csv', n_shuffle=1000, n_rewire=100,
                    proteome_id='UP000005640'): # refer to runNEE
    prot = pd.read_csv(proteome_file, sep='\t', keep_default_na=False)
    prot = prot[prot['Proteomes'].str.contains(proteome_id)].reset_index(drop=True) # same entries as NEE_decomposition.py
    indx = {e: i for i, e in enumerate(prot.Entry)}
    ints = read_interactions(pickle_file)
    pairs = {(indx[a], indx[b]) for a in ints for b in ints[a] if a in indx and b in indx} # pooled partner links (both directions)
    loops = sorted(i for i, j in pairs if i == j) # self-interactions
    edges = sorted({(min(i, j), max(i, j)) for i, j in pairs if i != j}) # undirected interactions between different proteins
    n = len(prot)
    def adjacency(edge_list): # symmetric partner matrix with the self-interactions
        r = [i for i, j in edge_list] + [j for i, j in edge_list] + loops
        c = [j for i, j in edge_list] + [i for i, j in edge_list] + loops
        return sparse.csr_matrix((np.ones(len(r)), (r, c)), shape = (n, n))
    adj = adjacency(edges)
    rep = residue_counts(prot.Sequence.tolist())
    length = rep.sum(axis = 1, keepdims = True)

    data = pd.read_csv(data_file)
    data = data[data.Category.notna()] # analysed proteins
    rows = data.Entry.map(indx).values # their rows in the proteome table
    y, cat = data.Disease.values, data.Category.values
    out = [{'Null model': 'observed', 'Replicate': 0, 'AUC': roc_auc_score(y, nee_of(rep, adj)[rows]),
            'tau_b': stats.kendalltau(nee_of(rep, adj)[rows], cat, variant='b')[0]}]

    rng = np.random.default_rng(2026) # random numbers of the shuffles
    for k in range(n_shuffle):
        perm = rng.permutation(n) # protein i gets the composition of protein perm[i]
        rep_shuffled = rep[perm] / length[perm] * length # scaled to the length of protein i
        nee = nee_of(rep_shuffled, adj)[rows]
        out.append({'Null model': 'composition shuffle', 'Replicate': k + 1, 'AUC': roc_auc_score(y, nee),
                    'tau_b': stats.kendalltau(nee, cat, variant='b')[0]})

    random.seed(2026) # igraph uses the random numbers of Python
    for k in range(n_rewire):
        graph = igraph.Graph(n=n, edges=edges)
        graph.rewire(n=10 * len(edges), allowed_edge_types='simple') # ten switch attempts per interaction
        nee = nee_of(rep, adjacency([tuple(e) for e in graph.get_edgelist()]))[rows]
        out.append({'Null model': 'degree-preserving rewiring', 'Replicate': k + 1, 'AUC': roc_auc_score(y, nee),
                    'tau_b': stats.kendalltau(nee, cat, variant='b')[0]})
    out = pd.DataFrame(out)
    out.to_csv(outfile, index=False) # writing the CSV output
    print(out.groupby('Null model')[['AUC', 'tau_b']].agg(['mean', 'std', 'min', 'max']).round(4).to_string())
    return out

