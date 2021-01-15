#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:49:31 2021

@author: babemomo
"""
import math as m


def ProfileMatrix(kmer_list, k, t, mode):  # nrow = t, ncol = k
    """Input a list of kmers and output a entropy score or profile probilities.

    It can be switched to hamming distance score following #.
    """
    matrix = [[0] * k for i in range(4)]  # initiate a zero-filled nested list as matrix.
    for col, colbases in enumerate(zip(*kmer_list)):
        for row, letter in enumerate('ACGT'):
            matrix[row][col] = colbases.count(letter)
    # hamming score
    if mode == 'ham':
        return t * k - sum([max(col) for col in zip(*matrix)])
    pseudo_count = 0.1
    for row in range(4):
        for col in range(k):
            matrix[row][col] += pseudo_count
            matrix[row][col] /= (t + 4 * pseudo_count)
    if mode == 'profile':
        return matrix
    if mode == 'score':
        # entropy score = sigma(-p*log_2(p)). Max score of 'ACGT' colume is 2.
        score = -sum([matrix[row][col] * m.log2(matrix[row][col])
                     for col in range(k) for row in range(4)])
        return score


rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def ProfileMostPKmer(sequence, k, matrix):
    """Input a sequence and a profile matrix (np_matrix), a score-optimized
    motif is returned.

    Loop through all substrings of the sequence.
    """
    max_p = 0
    motif = []
    for i in range(len(sequence) - k + 1):
        p = 1
        kmer = sequence[i: (i + k)]
        for col, letter in enumerate(kmer):
            p *= matrix[rows[letter]][col]
        if p > max_p:
            max_p = p
            motif = [kmer]
        elif p == max_p:
            motif.append(kmer)
    return motif[0]
