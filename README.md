# write-counts

Version 0.1.0

Python script to write a single gene count table from multiple count files. Repo includes compressed example data and reference files.

## Overview

Input:
- One .hdf5 file per sample (can be many samples, grouped by sequencing LIMSID, or other variable of choice).

The filenames *must* contain unique combinations of identifiers, that are also found as column names in the metadata file. The hdf5 files contain vectors of counts, per gene ID (with gene IDs sorted lexicographically).

Output:
- One table (genes x samples) per group.

Requires metadata:
- One table (samples x variables) per group. 

The metadata table must contain columns with names that correspond with the filenames (because files will be loaded based on the metadata annotations).

Requires reference genome file to extract gene IDs:
- Homo_sapiens.GRCh37.dna.gencode.v26lift37.basic.with_ERCC.gff 

This file is too large to host on GitHub without compression, so be sure to unarchive after cloning the repo.

### Workflow
- Clone the repo from main
- Verify folder structure
- ! Unarchive /ref gff file
- ! Unarchive /data/raw/ folder
- Navigate to root (where the main Python script is located)
- Verify and install dependencies (in particular gtfparser; also see requirements.txt)
- Run script

## Project organization

```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── write_count_table.py      <- Main Python script.
├── data                      <- 
│   ├── processed             <- The final count table.
│   ├── raw                   <- Count files in hdf5 format. Decompress before use.
├── docs                      <- Documentation notebook for users 
│   ├── manuscript            <- Manuscript source, e.g., LaTeX, Markdown, etc. 
│   └── reports               <- Other project reports and notebooks (e.g. Jupyter, .Rmd)
├── metadata                  <- Any annotations to the data. 
├── ref                       <- Reference gene annotation file. Decompress before use.
├── results
│   ├── figures               <- Figures for the manuscript or reports (PG)
│   └── output                <- Other output for the manuscript or reports 
└── src                       <- Source code for this project 

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

## Citation

Please [cite this project as described here](/CITATION.md).
