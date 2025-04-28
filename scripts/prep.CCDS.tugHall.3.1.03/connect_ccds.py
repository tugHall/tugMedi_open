#!/usr/bin/env python3
# -*- coding: utf-8; tab-width: 4 -*-


import re
import sys
from collections import defaultdict
from pprint import pprint

import pandas as pd
from icecream import ic

ic.configureOutput(includeContext=True)

ic.disable()

def convert_locs(locs_str):
    ic(locs_str)
    locs_str = re.sub(r'[\s\[\]]', '', locs_str)
    locs_str = re.sub(',$', '', locs_str)
    locs = []
    for entry in re.split(',', locs_str):
        entry = entry.strip()
        start, end = re.split('-', entry)
        if re.match(r'^[\d]+$', start) and re.match(r'^[\d]+$', end):
            locs.append({'start': int(start), 'end': int(end)})
    return(locs)


def get_pseudo_gene(chrom, start, end, locs):
    cds_location = ', '.join(['{}-{}'.format(l['start'], l['end']) for l in locs])
    cds_location = '[' + cds_location + ']'

    return(pd.DataFrame({
        '#chromosome': [chrom],
        'nc_accession': ['NC_allGenes_chr' + chrom],
        'gene': ['_allGenes_chr' + chrom],
        'gene_id': [chrom],
        'ccds_id': ['CCDS_allGenes_chr' + chrom],
        'ccds_status': ['Public'],
        'cds_strand': ['+'],
        'cds_from': [str(start)],
        'cds_to': [str(end)],
        'cds_locations': [cds_location],
        'match_type': ['Identical'],
        'cds_len': [str(end - start)]
    }))
   

def check_overlap(locs):
    locs = sorted(locs, key=lambda x:(x['start'], x['end']))
    uniq_locs = []
    prev = {'start': 0, 'end': 0}
    start = 0
    end = 0
    for loc in locs:
        if loc['start'] <= prev['start'] <= loc['end'] \
           or loc['start'] <= prev['end'] <= loc['end']:
            pass
        else:
            uniq_locs.append(loc)
            if start == 0 or start > loc['start']:
                start = loc['start']
            if end == 0 or end < loc['end']:
                end = loc['end']
        prev = loc
    return(start, end, uniq_locs)

### main

ccds_input_file = sys.argv[1]
ccds_output_file = sys.argv[2]
ccds_gene_file = sys.argv[3]

gene_df = pd.read_csv(ccds_gene_file, header=0, index_col=False, sep='\t', 
                      keep_default_na=False, na_values=None, dtype='object')
genes = gene_df['NCBI gene ID']

df = pd.read_csv(ccds_input_file, header=0, index_col=False, sep='\t', 
                 keep_default_na=False, na_values=None, dtype='object')

df = df[df['ccds_status'] == 'Public']

df_target = df[df['gene_id'].isin(genes)].copy()
df_other = df[~df['gene_id'].isin(genes)].copy()
del df

### target genes

df_target[['ccds_id.major', 'ccds_id.suf']] = df_target.apply(
    lambda x: 
    pd.Series(x['ccds_id'].split('.', 1))
    if '.' in x['ccds_id']
    else pd.Series([x['ccds_id'], '0']), 
    axis=1
)

# CCDS ex. CCDS15.1
df_target['ccds_id.major'] = df_target.apply(
    lambda x: re.sub(r'^CCDS', r'', x['ccds_id.major']), 
    axis=1
)

df_target = df_target.astype(
    {'ccds_id.major': int, 'ccds_id.suf': int }, 
    errors='raise'
)

df_target.sort_values(['gene_id', 'ccds_id.major', 'ccds_id.suf'], 
                      inplace=True)
df_target.drop_duplicates(subset=['gene_id'], keep='first', inplace=True)

### pseudo genes

for col in ['cds_to', 'cds_from']:
    df_other[col] = df_other[col].str.replace('-', '0').astype(int)
df_other = df_other.astype({'cds_to': int, 'cds_from': int}, errors='raise')

df_other['cds_len'] = df_other['cds_to'] - df_other['cds_from']
df_group = df_other.groupby('gene')
df_max = df_other.loc[df_group['cds_len'].idxmax(), :]

df_list = []
for chrom, group in df_other.groupby('#chromosome'):
    if chrom.isdecimal(): ## execpt for chrom X, Y
        locs = []
        for index, row in group.iterrows():
            locs.extend(convert_locs(row['cds_locations']))

        start, end, uniq_locs = check_overlap(locs)
        df_list.append(get_pseudo_gene(chrom, start, end, uniq_locs))

df_output = pd.concat(df_list)
df_output.drop(columns=['cds_len'], inplace=True)

### merge pseudo genes + target genes

df_output = pd.concat([df_output, df_target[df_output.columns]], axis=0)
df_output.to_csv(ccds_output_file, sep='\t', index=False)
