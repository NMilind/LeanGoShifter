# Lean GoShifter

[GoShifter](https://github.com/immunogenomics/goshifter) is an analysis tool to test for enrichment of GWAS loci in genomic annotations. For details, please refer to the paper.

Trynka G, Westra H-J, Slowikowski K, Hu X, Xu H, Stranger BE, et al. Disentangling the Effects of Colocalizing Genomic Annotations to Functionally Prioritize Non-coding Variants within Complex-Trait Loci. The American Journal of Human Genetics. 2015 Jul 2;97(1):139–52. 
[https://doi.org/10.1016/j.ajhg.2015.05.016](https://doi.org/10.1016/j.ajhg.2015.05.016).

I have reimplemented the tool in Python 3 with a reduced feature set:

1. Rather than calculating SNPs in LD with the lead SNP at a locus, Lean GoShifter expect the user to have pre-calculated a fine-mapped SNP set for each locus.
2. Lean GoShifter does not implement the conditional analysis for colocalizing annotations.

## Requirements

Lean GoShifter is written for **Python 3.8**. It uses the following base modules:

1. `argparse`
2. `functools`
3. `math`
4. `multiprocessing`
5. `os`
6. `sys`

In addition, Lean GoShifter requires the following modules:

1. `numpy >= 1.21.3`
2. `pandas >= 1.3.4`

You can install these using the following command:

```
python3 -m pip install numpy==1.21.3 pandas==1.3.4
```

## How to use Lean GoShifter

### Description

The script takes as input a list of SNPs assigned to each locus. These may be the lead SNP and tagging LD SNPs or a fine-mapped set of SNPs around the lead SNP. The script also takes a set of genomic intervals that represent annotations.

Loci are calculated as the furthest SNPs at each locus extended by twice the median size of the annotation. The annotations at each locus are permuted and the enrichment of overlapping SNPs at loci is estimated using the null distribution. Refer to the original paper for specific details.

### Usage

```
usage: lean_go_shifter.py [-h] [--threads THREADS] snp_map annotation permute out_dir prefix

A lean reimplementation of GoShifter

positional arguments:
  snp_map            A table of SNPs for each individual signal.
  annotation         File with genomic regions representing annotations.
  permute            Number of permutations for the null distribution.
  out_dir            Output directory.
  prefix             Prefix for output files.

optional arguments:
  -h, --help         show this help message and exit
  --threads THREADS  The number of threads to use.
```

### Inputs

This is the expected format for the `snp_map` argument. The file is tab-delimited. The `Locus` column indicates separate signals to test for enrichment in the genomic annotations. This example contains 2 loci with a fine-mapped set of one SNP and a third locus with a fine-mapped set of 7 SNPs.

```
Locus	SNP	Chr	Position
1	rs3131972	1	817341
2	rs2272757	1	946247
3	rs13302945	1	953779
3	rs13302957	1	955641
3	rs13302996	1	966227
3	rs13303051	1	954333
3	rs13303056	1	953778
3	rs13303206	1	954258
3	rs13303227	1	957365
```

This is the expected format for the `annotation` argument. The file is tab-delimited.

```
Chr	Start	End
1	23224	23849
1	24025	26778
1	28524	30136
1	30354	30875
1	31109	31894
1	32151	32330
1	32460	33367
1	33645	33887
1	33988	35232
```

### Example

```
mkdir test_output/
python3 lean_go_shifter.py snp_map.tsv annotation.tsv 1000 test_output/ prefix
```

### Outputs

1. `test_output/prefix_p_value.tsv` - Contains the p-value of enrichment of loci overlapping the genomic annotation.
2. `test_output/prefix_overlap_scores.tsv` - Contains the loci that were tested for enrichment. A column indicates if the locus overlapped any genomic annotations. Another column indicates the likelihood of an overlap occuring, with low scores representing unlikely overlap that drives the enrichment of the genomic annotation. Refer to the paper for specific details of the overlap score.
