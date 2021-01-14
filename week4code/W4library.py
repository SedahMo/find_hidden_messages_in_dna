#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 00:53:00 2021

@author: babemomo
"""
import numpy as np


rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def ProfileMatrix(kmer_list, k, t, mode):  # nrow = t, ncol = k
    """Input a list of kmers and output a entropy score or profile probilities.

    It can be switch to hamming distance score following #.
    """
    np_matrix = np.zeros((4, k))  # matrix.shape = 4, k
    for col, colbases in enumerate(zip(*kmer_list)):
        for letter in 'ACGT':
            np_matrix[rows[letter], col] = colbases.count(letter)
    pseudo_count = 0.1
    if mode == 'profile':
        np_matrix += 1 * pseudo_count
        return np_matrix / (t + 4 * pseudo_count)
    # if mode == 'score':
    #     return t * k - sum((np.amax(np_matrix, axis=0)))
    if mode == 'score':
        np_matrix += 1 * pseudo_count
        np_matrix /= (t + 4 * pseudo_count)
        # entropy score = sigma(-plog_2(p)). Max score of 'ACGT' colume is 2.
        return np.sum(-np_matrix * np.log2(np_matrix))
    if mode == 'entropy':
        np_matrix /= t
        np_matrix = np_matrix.reshape(-1)
        np_matrix = np_matrix[np_matrix != 0]
        return np.sum(-np_matrix * np.log2(np_matrix))


def ProfileMostPKmer(sequence, k, np_matrix):
    """Input a sequence and a profile matrix (np_matrix), a score-optimized
    motif is returned.

    Loop through all substrings of the sequence."""
    max_p = 0
    motif = []
    for i in range(len(sequence) - k + 1):
        p = 1
        kmer = sequence[i: (i + k)]
        for col, letter in enumerate(kmer):
            p *= np_matrix[rows[letter], col]
        if p > max_p:
            max_p = p
            motif = [kmer]
        elif p == max_p:
            motif.append(kmer)
    return motif[0]
