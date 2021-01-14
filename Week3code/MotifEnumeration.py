#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:10:41 2021

@author: Monica Ruan
"""

import W3library as w3


def MotifEnumeration(dna_list, k, d):
    pattern_set = None
    for sequence in dna_list:
        # Generate all the neighbors of sequences in dna_list.
        neighbor_set_temp = set()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i: (i + k)]
            neighbor_set_temp |= w3.IterativeNeighbors(kmer, d)
        # Calculate the intersection of all neighbors_set.
        if pattern_set is not None:
            pattern_set &= neighbor_set_temp
        else:
            pattern_set = neighbor_set_temp
    return pattern_set


if __name__ == '__main__':
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        params = line_list.pop(0)
        dna_list = line_list
        k = int(params[0])
        d = int(params[2])
    print(*MotifEnumeration(dna_list, k, d))
