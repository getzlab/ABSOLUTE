# ABSOLUTE

Estimate tumor purity, ploidy, and absolute copy number from a copy number segmentation.

Publication: [Absolute quantification of somatic DNA alterations in human cancer (Carter et al. 2012)](https://doi.org/10.1038/nbt.2203).

[This page on the GenePattern website contains some helpful information as well.](https://www.genepattern.org/analyzing-absolute-data#gsc.tab=0)

## High-level concepts

ABSOLUTE receives the following inputs:

* Segmented, smoothed (relative) copy number profiles
    - I.e., from HapASeg or AllelicCapSeg
* Somatic point mutations with their allelic fractions

ABSOLUTE uses these as evidence to infer the following quantities:

* Tumor purity
* Tumor ploidy
* Absolute copy numbers for genome segments

## Outputs from ABSOLUTE


The ABSOLUTE R script (`ABSOLUTE_cli_start.R`) outputs the following files:

* `*.ABSOLUTE_plot.pdf`. Plots of ABSOLUTE candidate solutions
* `*.PP-modes.data.RData`. RData file containing candidate solutions

Once a candidate solution is chosen, the "extractor" script (`ABSOLUTE_extract_cli_start.R`) can produce the following outputs:
* `{sample_name}_ABS_MAF.txt`
* `{sample_name}.segtab.txt`
* `{sample_name}.ABSOLUTE.{analyst_id}.called.RData`
    - RData file containing a list of candidate ABSOLUTE solutions
* `{sample_name}.{analyst_id}.ABSOLUTE.table.txt`
    - Contains purities and ploidies. Analyst can edit this file with manual annotations in order to

