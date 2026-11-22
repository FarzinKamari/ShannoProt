#!/usr/bin/env python
# coding: utf-8

# # GTDB genomes of the proteomes -- phylogenetic signal and genome quality (rev 26)
# * Each bacterial and archaeal proteome is matched to a genome of the Genome Taxonomy Database (GTDB, release 232): through the NCBI assembly that gave its GC content (GC null per proteome.csv; the same accession, another version of it, or the GenBank accession of a RefSeq genome), otherwise through its NCBI taxid (the most frequent GTDB species representative among the GTDB genomes of the taxid).
# * One row per GTDB species representative, with the median CEI, GC, expected CEI and delta CEI of its proteomes: the input of "Phylogenetic signal.R" (Pagel's lambda and PGLS on the GTDB trees).
# * The proteomes above the highest CEI of the uniform-GC null, by GC assignment, assembly level, RefSeq category, genome quality (CheckM2 completeness and contamination in the GTDB metadata) and genome origin (metagenome-assembled, environmental or single-cell genomes), and whether they lie below the codon-position null of their genus.
# ##### Input files:
# >* GC null per proteome.csv and GC null maxima.csv (from Proteome_GC_null.py), Codon position GC per genus.csv (from CDS_codon_position_GC.py; optional)
# >* bac120_metadata_r232.tsv.gz and ar53_metadata_r232.tsv.gz (GTDB release 232)
# ##### Output files:
# >* Proteomes GTDB.csv (one row per proteome), GTDB species values.csv (one row per GTDB species) and GC ceiling by genome quality.csv

# In[1]:


import os
import numpy as np
import pandas as pd
from CDS_codon_position_GC import genus_of


# In[2]:


def gtdb_matching(gc_file='GC null per proteome.csv', metadata_files=('bac120_metadata_r232.tsv.gz', 'ar53_metadata_r232.tsv.gz'),
                  outfile='Proteomes GTDB.csv'): # refer to runGC
    cols = ['accession', 'ncbi_genbank_assembly_accession', 'gtdb_genome_representative', 'gtdb_taxonomy', 'checkm2_completeness',
            'checkm2_contamination', 'ncbi_genome_category', 'ncbi_taxid', 'ncbi_assembly_level']
    meta = pd.concat([pd.read_csv(f, sep='\t', usecols=lambda c: c in cols, low_memory=False) for f in metadata_files], ignore_index=True)
    meta['acc'] = meta.accession.str[3:] # RS_GCF_... and GB_GCA_... -> GCF_... and GCA_... (with version)
    meta['base'] = meta.acc.str.split('.').str[0] # without version
    meta['gb_base'] = meta.ncbi_genbank_assembly_accession.astype(str).str.split('.').str[0] # GenBank accession of a RefSeq genome
    gc = pd.read_csv(gc_file, low_memory=False)
    gc['base'] = gc['Assembly accession'].astype(str).str.split('.').str[0]
    lookup = {col: meta.drop_duplicates(col).set_index(col) for col in ('acc', 'base', 'gb_base')}
    rows = []
    for _, r in gc.iterrows():
        hit, how = None, ''
        acc = r['Assembly accession']
        if not isinstance(acc, str): # no assembly with a GC value
            rows.append((None, ''))
            continue
        if acc in lookup['acc'].index:
            hit, how = lookup['acc'].loc[acc], 'accession'
        elif r.base in lookup['base'].index:
            hit, how = lookup['base'].loc[r.base], 'accession (other version)'
        elif r.base.startswith('GCA') and r.base in lookup['gb_base'].index:
            hit, how = lookup['gb_base'].loc[r.base], 'GenBank accession of a RefSeq genome'
        rows.append((hit['accession'] if hit is not None else None, how))
    gc['GTDB genome'], gc['GTDB match'] = zip(*rows)
    # otherwise by taxid: the most frequent species representative among the GTDB genomes of the taxid
    tx = meta.groupby('ncbi_taxid').gtdb_genome_representative.agg(lambda s: s.value_counts().index[0])
    miss = gc['GTDB genome'].isna() & gc.Domain.isin(['Bacteria', 'Archaea'])
    taxid = gc['Current taxid'].fillna(gc.Organism_ID).astype(int)
    gc.loc[miss, 'GTDB species representative'] = taxid[miss].map(tx)
    gc.loc[miss & gc['GTDB species representative'].notna(), 'GTDB match'] = 'taxid'
    m = meta.set_index('accession')
    found = gc['GTDB genome'].notna()
    gc.loc[found, 'GTDB species representative'] = gc.loc[found, 'GTDB genome'].map(m.gtdb_genome_representative)
    for c, name in (('checkm2_completeness', 'CheckM2 completeness'), ('checkm2_contamination', 'CheckM2 contamination'),
                    ('ncbi_genome_category', 'NCBI genome category')):
        gc[name] = gc['GTDB genome'].map(m[c]) # genome quality and origin: genomes matched by accession only
    gc['GTDB taxonomy'] = gc['GTDB species representative'].map(m.gtdb_taxonomy)
    gc = gc.drop(columns='base')
    gc.to_csv(outfile, index=False) # writing the CSV output
    pro = gc[gc.Domain.isin(['Bacteria', 'Archaea'])]
    print(len(pro), 'prokaryotic proteomes matched to GTDB:', pro['GTDB match'].replace('', np.nan).fillna('not matched').value_counts().to_dict())
    return gc


# In[3]:


def gtdb_species_values(mapping_file='Proteomes GTDB.csv', outfile='GTDB species values.csv'): # refer to runGC
    g = pd.read_csv(mapping_file, low_memory=False)
    g = g[g['GTDB species representative'].notna() & g.GC.notna()]
    tips = g.groupby('GTDB species representative').agg(Domain=('Domain', 'first'), Proteomes=('CEI', 'size'), CEI=('CEI', 'median'),
                                                         GC=('GC', 'median'), CEI_null=('CEI null', 'median'), Delta_CEI=('Delta CEI', 'median'),
                                                         Taxonomy=('GTDB taxonomy', 'first')).reset_index() # median of the proteomes of each species
    tips = tips.rename(columns={'CEI_null': 'CEI null', 'Delta_CEI': 'Delta CEI', 'Taxonomy': 'GTDB taxonomy'})
    tips.to_csv(outfile, index=False) # writing the CSV output
    print(len(tips), 'GTDB species:', tips.Domain.value_counts().to_dict(), '| proteomes per species: median', tips.Proteomes.median(),
          'max', tips.Proteomes.max())
    return tips


# In[4]:


def gc_ceiling_quality(mapping_file='Proteomes GTDB.csv', maxima_file='GC null maxima.csv', codon_file='Codon position GC per genus.csv',
                       outfile='GC ceiling by genome quality.csv'): # refer to runGC
    g = pd.read_csv(mapping_file, low_memory=False)
    ceiling = float(pd.read_csv(maxima_file).iloc[0]['Highest CEI']) # the highest CEI of the uniform-GC null (genetic code 11)
    g['above'] = g.CEI > ceiling
    g['high quality'] = (g['CheckM2 completeness'] >= 95) & (g['CheckM2 contamination'] <= 5)
    g['MAG'] = g['NCBI genome category'].astype(str).str.contains('metagenome|environmental|single cell', case=False)
    rows = []

    def share(x, domain, grouping, category):
        rows.append({'Domain': domain, 'Grouping': grouping, 'Category': category, 'Above the maximum': int(x.above.sum()),
                     'Proteomes': len(x), 'Percent above': round(100 * x.above.mean(), 1)})

    for dom in ('Bacteria', 'Eukaryota', 'Archaea'):
        x = g[(g.Domain == dom) & g.GC.notna()] # proteomes with a genome GC
        share(x, dom, 'All proteomes with a genome GC', 'all')
        for c in ('Match type', 'Assembly level', 'RefSeq category'):
            for k, v in x.groupby(c):
                share(v, dom, c, k)
        if dom != 'Eukaryota': # GTDB holds bacteria and archaea only
            gx = x[x['GTDB genome'].notna()]
            for k, v in gx.groupby('high quality'):
                share(v, dom, 'GTDB genome quality', 'CheckM2 >= 95% complete and <= 5% contamination' if k else 'other')
            for k, v in gx.groupby('MAG'):
                share(v, dom, 'Genome origin', 'metagenome-assembled, environmental or single cell' if k else 'other')
            best = x[(x['Match type'] == 'taxid') & x['Assembly level'].isin(['Complete Genome', 'Chromosome']) & x['high quality'] & ~x['MAG']]
            share(best, dom, 'Strictest subset', 'GC from the same taxid, complete or chromosome-level, high quality, not metagenome-assembled')
    out = pd.DataFrame(rows)
    if os.path.exists(codon_file): # codon-position GC of the genus (one annotated genome per genus): below the position-specific null?
        cp = pd.read_csv(codon_file)
        g['Genus'] = g.Organism.apply(genus_of)
        pos = g.merge(cp[['Genus', 'CEI null (GC1, GC2, GC3)']], on='Genus', how='left')
        for dom in ('Bacteria', 'Eukaryota', 'Archaea'):
            x = pos[(pos.Domain == dom) & pos.above & pos['CEI null (GC1, GC2, GC3)'].notna()] # above the maximum, genus with a null
            out = pd.concat([out, pd.DataFrame([{'Domain': dom, 'Grouping': 'Codon-position null of the genus',
                                                 'Category': 'above the maximum, genus with a codon-position null', 'Above the maximum': len(x),
                                                 'Proteomes': len(x), 'Percent above': 100.0,
                                                 'Below the codon-position null': int((x.CEI <= x['CEI null (GC1, GC2, GC3)']).sum())}])],
                            ignore_index=True)
    out.to_csv(outfile, index=False) # writing the CSV output
    print(out.to_string(index=False))
    return out
