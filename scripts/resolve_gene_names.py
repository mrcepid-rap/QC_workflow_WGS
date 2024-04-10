#!/usr/bin/env python3
import sys
import pandas as pd
import requests
import functools
import operator

from pathlib import Path

from requests import JSONDecodeError

"""This script takes the Path to an annotation .tsv file and reads it, filters to the requested chromosome, and
adds a gene symbol columns to allow for easier merging of alphafold annotations.

The use of gene symbols was chosen because alpha_missense uses UniProt or ENST IDs to denote the annotated
transcript. The ENST IDs often line-up with the MANE / Canonical transcript, but can also be IDs of unknown
origin; often these missing IDs are old ENSEMBL IDs that are no longer used and are difficult to process.

I then query the ENSEMBL REST API (https://rest.ensembl.org/) to get the ENSEMBL-chosen gene name for that
transcript which _should_ line up with what was used when running VEP to do the original variant annotation.

:param annotation_path: Path to a gene-specific annotation file
"""

annotation_path = Path(sys.argv[1])
annotation_path = Path('eve/eve.tsv.gz')

# Load the raw table into Pandas
annotation_table = pd.read_csv(annotation_path, sep='\t')

annotation_table.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

# Need a non-redundant list of genes to query
enst_values = list(set(annotation_table['ENST'].to_list()))
matched_frame = pd.DataFrame(data={'ENST': enst_values})
matched_frame['ENST_split'] = matched_frame['ENST'].apply(lambda x: set(x.split('|')))

gene_list = [list(x) for x in matched_frame['ENST_split'].to_list()]
gene_list = functools.reduce(operator.iconcat, gene_list, [])
gene_list = list(set(gene_list))

# Set current REST API information
server = 'http://rest.ensembl.org'
ext = "/lookup/id"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

archive_genes = {}
failed_genes = []
finalised_genes = {}

# REST API POST can only take 1k queries at a time, so have to split and run in chunks (I do smaller to make sure...)
for i in range(0, len(gene_list), 100):

    curr_gene_list = gene_list[i:i+100]

    gene_query = '{ "ids": ["' + '","'.join(curr_gene_list) + '"] }'
    curr_req = requests.post(server+ext, headers=headers, data=gene_query)
    gene_json = curr_req.json()

    for gene in curr_gene_list:

        curr_gene_json = gene_json[gene]

        if curr_gene_json is None:
            # archive_genes[gene.ENST] = gene.uniprot_id
            failed_genes.append(gene)
        else:
            if 'display_name' in curr_gene_json:
                finalised_genes[gene] = curr_gene_json['display_name'].split('-')[0]
            else:
                failed_genes.append(gene)

gene_table = pd.DataFrame(index=list(finalised_genes.keys()), data=finalised_genes.values(), columns=['SYMBOL'])
alpha_fold_gene_id = annotation_table.merge(gene_table, left_on='ENST', right_index=True, how='left')
