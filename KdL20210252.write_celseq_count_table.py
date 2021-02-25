#!/usr/bin/env python3

# Kim de Luca
# 20210219

#### Write count table from binned count (hdf5) files ####

# Load modules
import os
from glob import glob # unused
from itertools import chain, groupby

import h5py
import pandas as pd
import numpy as np
import sys # unused

# %pip install git+https://github.com/krooijers/gtfparser --user
from gtfparser import GTFParser

# Define LIMSID of interest
LIMSID = 'KIN3801'

# Define metadata files
ANNOFN = "./metadata/{LIMSID}_anno.tsv".format(LIMSID=LIMSID)
GTFFILE = "./metadata/Homo_sapiens.GRCh37.dna.gencode.v26lift37.basic.with_ERCC.gff"

# Define filename formats
#TODO fix manual edit of LIMSID below
FNFMT = "./data/raw/KIN3801.index{indexnr:02d}.BCv3set2_BC_{bcnr:03d}.counts.hdf5"
OUTFNFMT = "./data/processed/tables/{limsid}.counts_celseq.tsv.gz"

# Define features and attributes of interest
FEATURE = 'exon'
ATTR_KEY = 'gene_id'
COLFMT = "{limsid}.index{indexnr:02d}.BCv3set2_BC_{bcnr:03d}"

# parse GTF file to get list of gene names
gtf = GTFParser(GTFFILE)
gtf = filter(lambda f: f.feature == FEATURE, gtf)

chrom_strand_key = lambda f: (f.iv.chrom, f.iv.strand)
attr_groupby_key = lambda f: f.attr.get(ATTR_KEY)

gtf = chain.from_iterable(map(lambda t: t[1], groupby(gtf, key=chrom_strand_key)))
geneids = set()
geneid2name = dict()
for geneid, fiter in groupby(gtf, key=attr_groupby_key):
    f = next(fiter)
    name = f.attr.get("gene_name")
    geneid2name[geneid] = name
    geneids.add(geneid)

geneids = sorted(geneids)

# Load metadata file
anno = pd.read_table(ANNOFN, sep='\t')

countdata = list()
for irow, row in anno.iterrows():
    fn = FNFMT.format(limsid=row["limsid"], indexnr=row["indexnr"], bcnr=row["barcodenr"]) # change according to samplesheet
    if not os.access(fn, os.R_OK):
        print("Warning: no data for %s" % fn)
        # this happens if tophat2 didnt run because the demultiplexed input FASTQ files were empty...
        d = np.zeros(len(geneids), dtype=int)
    else:
        f = h5py.File(fn, 'r')
        d = f['counts'][:]
    countdata.append(d)

countdf = pd.DataFrame(np.array(countdata).T)
cols = [COLFMT.format(limsid=row["limsid"], indexnr=row["indexnr"], bcnr=row["barcodenr"]) for row in [anno.iloc[irow] for irow in range(len(countdata))]] # change according to samplesheet
countdf.columns = cols
countdf["gene_id"] = geneids
countdf["gene_name"] = [geneid2name[geneid] for geneid in geneids]
countdf.set_index("gene_id", drop=True, inplace=True)
countdf = countdf[["gene_name"] + sorted(cols)]

outfn = OUTFNFMT.format(limsid=row["limsid"])

# Write the count table to disk
countdf.to_csv(outfn, sep="\t", compression="gzip")
