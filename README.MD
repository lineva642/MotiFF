# MotiFF

An open source software for recognition of consensus motifs in proteomics
## Overview

**MotiFF** is a tool for large scale proteomics that discovers the enriched amino acid patterns containing modified amino acids in the results of proteomic analysis. MotiFF implements calculation of the binomial probabilities and χ2  criterion to identify statistically significant amino acid patterns or motifs for the modified peptides. The algorithm was implemented in Python as an open-source command-line tool available at https://github.com/lineva642/PTM.
## Installation

Download from Github repository https://github.com/lineva642/PTM. In the directory containing setup.py file run the following command:

    pip install .

## Usage

    MotiFF [dataset] [fasta] [-h] [--name_sample NAME_SAMPLE] [--interval_length INTERVAL_LENGTH] [--modification MODIFICATION] [--modification_site MODIFICATION_SITE] [--working_dir WORKING_DIR] [--algorithm ALGORITHM] [--p_value P_VALUE] [--occurrences OCCURRENCES] [-v {0,1,2}]
    
    positional arguments:
    dataset               Path to experimental dataset
    fasta                 Path to FASTA

    optional arguments:
    -h, --help            show this help message and exit
    --name_sample NAME_SAMPLE
                          Name of examined sample
    --interval_length INTERVAL_LENGTH
                          Number of amino acids before & after modified amino acid (default=6)
    --modification MODIFICATION
                          Name of modification(ex.PHOSPHORYLATION)
    --modification_site MODIFICATION_SITE
                          Modified amino acid (ex.S,T)
    --working_dir WORKING_DIR
                          Working dir for program (default=".")
    --algorithm ALGORITHM
                          Enter algorithm name: binom or chi2(default="chi2")
    --p_value P_VALUE     Enter p_value(default=0.000005)
    --occurrences OCCURRENCES
                          Enter number of motif occurrences in experimental dataset (default=20)
    -v {0,1,2}, --verbosity {0,1,2}
                          Output verbosity
                          
## Examples

## Output files
MotiFF produces summary table (motifs.csv) with discovered amino acid motifs and additional information:

| Column name | Description |
|:----------------|:----------------|
|motif|Enriched amino acid pattern containing modified amino acid|
|p_value|Probability of the identified motif|
|fg_matches|Number of modified peptides containing this enriched amino acid pattern|
|fg_size|Number of modified peptides among which search is carried out|
|bg_size|Number of FASTA peptides among which search is carried out|
|bg_matches|Number of FASTA peptides containing this enriched amino acid pattern|
