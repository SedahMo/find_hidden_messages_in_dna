#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 00:49:10 2021

@author: babemomo
"""

import W3library as w3
def HammingDistanceDiffLen(pattern, sequence):
    """Calculate the minimum Hamming distance between the pattern and every \
    kmer in a DNA sequence."""
    k = len(pattern)
    hamdistance = None
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i: (i + k)]
        if hamdistance is not None:
            hamdistance = min(hamdistance, w3.HammingDistance(pattern, kmer))
        else:
            hamdistance = w3.HammingDistance(pattern, kmer)
    return hamdistance


def MinHamDistance(pattern, dna_list):
    """Calculate the minimum Hamming distance from a DNA list."""
    return sum(HammingDistanceDiffLen(pattern, sequence) for sequence in
               dna_list)


if __name__ == "__main__":
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        pattern = line_list.pop(0)
        dna_list = line_list.pop(0).split()
    print(MinHamDistance(pattern, dna_list))
