#!/usr/bin/env python

'''
:File: lean_go_shifter.py
:Author: Nikhil Milind, Wellcome Sanger Institute, <nm18@sanger.ac.uk>
:Last Updated: 25 March 2022


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

BY USING THE SOFTWARE YOU ACKNOWLEDGE THAT YOU HAVE READ AND UNDERSTAND THE
TERMS OF USE BELOW. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

THIS SOFTWARE IS TO BE USED AS A RESEARCH TOOL ONLY. THE SOFTWARE TOOL SHALL
NOT BE USED AS A DIAGNOSTIC DECISION MAKING SYSTEM AND MUST NOT BE USED TO
MAKE A CLINICAL DIAGNOSIS OR REPLACE OR OVERRULE A LICENSED HEALTH CARE
PROFESSIONAL'S JUDGMENT OR CLINICAL DIAGNOSIS. ANY ACTION THE RESEARCHER TAKES
IN RESPONSE TO THE INFORMATION CONTAINED WITHIN IS AT THE RESEARCHER'S
DISCRETION ONLY.
'''


import argparse
import math
import os
import sys

import numpy as np
import pandas as pd


def parse_arguments():

    '''
    Parse and validate command line arguments.

    :return: A namespace containing the command line argument values.
    '''

    parser = argparse.ArgumentParser(description='A lean reimplementation of GoShifter')

    parser.add_argument('snp_map', help='A table of SNPs for each individual signal.')
    parser.add_argument('annotation', help='File with genomic regions representing annotations.')
    parser.add_argument('permute', help='Number of permutations for the null distribution.')
    parser.add_argument('out_dir', help='Output directory.')
    parser.add_argument('prefix', help='Prefix for output files.')

    args = parser.parse_args()

    if not os.path.isfile(args.snp_map):
        print(f'{args.snp_map} is not a valid path to a file', file=sys.stderr)
        exit(2)
    
    if not os.path.isfile(args.annotation):
        print(f'{args.annotation} is not a valid path to a file', file=sys.stderr)
        exit(2)
    
    try:
        args.permute = int(args.permute)
    except ValueError:
        print(f'{args.permute} is not an integer', file=sys.stderr)
        exit(2)
    
    if args.permute <= 0:
        print('Permutations must be greater than 0', file=sys.stderr)
        exit(2)

    return args


def read_snp_map(file_path):

    '''
    Reads the list of fine-mapped SNP sets from the flat file provided.

    :param file_path: File path to the SNP map.
    :return: A data frame of the snp map.
    '''

    snp_map = pd.read_csv(file_path, sep='\t', dtype={'Chr': str, 'Position': int})
    snp_map.sort_values(['Chr', 'Position'], ignore_index=True, inplace=True)

    return snp_map


def read_annotations(file_path):

    '''
    Reads the list of annotations from the flat file provided.

    :param file_path: File path to the annotations.
    :return: A data frame of the annotations.
    '''

    annotations = pd.read_csv(file_path, sep='\t', dtype={'Chr': str, 'Start': int, 'End': int})
    annotations.sort_values(['Chr', 'Start', 'End'], ignore_index=True, inplace=True)

    return annotations


def infer_loci(snp_map, annotations):

    '''
    Infers the genomic range for each locus. Takes the minimum and maximum positions provided
    per set of SNPs within a locus and adds the median width of an annotation to the ends.

    :param snp_map: Set of SNPs at each locus.
    :param annotations: Set of annotations over genome.
    :return: A data frame of ranges for the loci.
    '''

    median_annotation_width = math.floor((annotations['End'] - annotations['Start']).median())

    loci = snp_map.groupby('Locus').agg({
        'Chr': 'first',
        'Position': [
            lambda x: x.min() - 2 * median_annotation_width,
            lambda x: x.max() + 2 * median_annotation_width
        ]
    })
    loci.columns = loci.columns.droplevel(0)
    loci.columns = ['Chr', 'Start', 'End']
    loci.reset_index(inplace=True)
    loci.sort_values(['Chr', 'Start', 'End'], ignore_index=True, inplace=True)

    return loci


def snps_overlap_with_annotations(snp_map, annotations):

    '''
    Tests if any SNP in the list overlaps any annotation.

    :param snp_map: The list of SNPs to test for overlap.
    :param annotations: The annotations to test for overlap.
    :return: True if any SNP overlaps any annotation.
    '''

    annotation_index = 0
    snp_index = 0

    annotation_iterator = annotations.iterrows()
    snp_iterator = snp_map.iterrows()

    annotation_index, annotation = next(annotation_iterator, (None, None))
    snp_index, snp = next(snp_iterator, (None, None))

    # Iterate until either the annotations or SNPs are all considered for overlap
    while annotation_index is not None and snp_index is not None:

        # If the annotation occurs before the SNP by genomic coordinates, increment the annotation
        # If the SNP occurs before the annotation by genomic coordinates, increment the SNP
        if annotation.Chr == snp.Chr:
            if annotation.Start <= snp.Position <= annotation.End:
                return True
                snp_index, snp = next(snp_iterator, (None, None))
            elif snp.Position < annotation.Start:
                snp_index, snp = next(snp_iterator, (None, None))
            else:
                annotation_index, annotation = next(annotation_iterator, (None, None))
        elif annotation.Chr < snp.Chr:
            annotation_index, annotation = next(annotation_iterator, (None, None))
        else:
            snp_index, snp = next(snp_iterator, (None, None))

    return False


def permute_locus(locus, snp_map, annotations, iters):

    '''
    Calculates the observed overlap with annotations and then performs permutations to
    compute a null distribution for the locus.

    :param locus: The bounds of the locus.
    :param snp_map: The list of SNPs at the locus.
    :param annotations: The annotations at the locus to permute.
    :param iters: The number of permutations to build the null distribution.
    :return: A tuple containing the observed overlap and the null distribution.
    '''

    if len(annotations) == 0:
        return 0, np.zeros((iters,))

    # Restrict annotations that are outside the locus bounds to strictly within the locus region
    annotations['Start'] = annotations['Start'].apply(lambda x: np.max((locus.Start, x)))
    annotations['End'] = annotations['End'].apply(lambda x: np.min((locus.End, x)))

    # Determine the observed overlap at the locus
    observed = int(snps_overlap_with_annotations(snp_map, annotations))

    # Build null distribution
    locus_width = locus.End - locus.Start + 1
    iterations = list()

    for i in range(iters):

        # Shift a random amount within the locus
        shift = np.random.randint(0, locus_width)

        # Add shift value to all genomic coordinates of the annotation
        annotations_iter = annotations.copy()
        annotations_iter[['Start', 'End']] = annotations_iter[['Start', 'End']].add(shift)

        # Identify any annotations that overflow (go past the boundary of the locus)
        annotations_iter_overflow = annotations_iter.copy()
        annotations_iter_overflow = annotations_iter_overflow.loc[annotations_iter_overflow['End'] > locus.End, :]

        # "Circularize" by mapping overflowing elements to the beginning of the locus
        annotations_iter_overflow['Start'] = locus.Start + (annotations_iter_overflow['Start'] - locus.End - 1)
        annotations_iter_overflow['End'] = locus.Start + (annotations_iter_overflow['End'] - locus.End - 1)

        # Trim overflowing elements to be within the locus
        annotations_iter['End'] = annotations_iter['End'].apply(lambda x: np.min((locus.End, x)))
        
        # Add back overflow segments and sort annotations for overlap test
        annotations_iter = pd.concat((annotations_iter_overflow, annotations_iter))
        annotations_iter.sort_values(['Chr', 'Start', 'End'], ignore_index=True, inplace=True)

        # Store iteration result
        iter_result = int(snps_overlap_with_annotations(snp_map, annotations_iter))
        iterations.append(iter_result)

    return observed, np.array(iterations)


def main():

    # Read command line arguments
    args = parse_arguments()

    # Read SNP map and annotation
    snp_map = read_snp_map(args.snp_map)

    annotations = read_annotations(args.annotation)

    # Infer loci from SNP map and annotations
    loci = infer_loci(snp_map, annotations)

    # Process each locus
    observed_dist = np.zeros((len(loci),))
    null_dist = np.zeros((args.permute, len(loci)))

    for locus_index, locus in loci.iterrows():

        # Identify SNPs and annotations at the locus
        snps_at_locus = snp_map.copy().loc[locus.Locus == snp_map['Locus'], :]
        snps_at_locus.reset_index(drop=True, inplace=True)

        annotations_at_locus = annotations.copy().loc[(locus.Chr == annotations['Chr']) & (locus.Start <= annotations['End']) & (annotations['Start'] <= locus.End), :]
        annotations_at_locus.reset_index(drop=True, inplace=True)
        
        # Perform permutations
        observed, iterations = permute_locus(locus, snps_at_locus, annotations_at_locus, args.permute)

        # Update observed and null distributions
        observed_dist[locus_index] = observed
        null_dist[:, locus_index] = iterations
    
    # Calculate statistics over count distribution
    count_observed = observed_dist.sum()

    count_null_dist = null_dist.sum(axis=1)

    p_value = (count_null_dist >= count_observed).sum() / args.permute

    with open(os.path.join(args.out_dir, f'{args.prefix}_p_value.tsv'), 'w') as f_out:
        f_out.write(f'pvalue\t{p_value}')

    # Calculate overlap scores for each locus
    overlap_scores = null_dist.sum(axis=0) / args.permute

    loci['Overlap'] = observed_dist
    loci['Overlap_Score'] = overlap_scores
    loci.sort_values(['Overlap', 'Overlap_Score'], ascending=[False, True], inplace=True)

    loci.to_csv(os.path.join(args.out_dir, f'{args.prefix}_overlap_scores.tsv'), sep='\t', index=False)


if __name__ == '__main__':

    main()
