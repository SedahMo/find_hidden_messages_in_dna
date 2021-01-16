#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:49:31 2021

@author: babemomo
"""
import math as m
import random as rd


def ProfileMatrix(kmer_list, k, t, mode, pseudo_count=1):  # nrow = t, ncol = k
    """Input a list of kmers and output a entropy score or profile probilities.

    It can be switched to hamming distance score by setting mode.
    """
    matrix = [[0] * k for i in range(4)]  # Initiate a zero-filled nested list.
    for col, colbases in enumerate(zip(*kmer_list)):
        for row, letter in enumerate('ACGT'):
            matrix[row][col] = colbases.count(letter)
    # hamming score
    if mode == 'ham':
        return t * k - sum([max(col) for col in zip(*matrix)])
    for row in range(4):
        for col in range(k):
            matrix[row][col] += pseudo_count
            matrix[row][col] /= (t + 4 * pseudo_count)
    if mode == 'profile':
        return matrix
    if mode == 'entr':
        # entropy score = sigma(-p*log_2(p)). Max score of 'ACGT' colume is 2.
        score = -sum([matrix[row][col] * m.log2(matrix[row][col])
                     for col in range(k) for row in range(4)])
        return score


rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def ProfileMostPKmer(sequence, k, matrix):
    """Input a sequence and a profile matrix, output a score-optimized motif.

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


def GibbsKmer(sequence, k, matrix):
    """Input a sequence and profile matrix, output a random choice of kmer.
    
    The kmer is based on gibbs sampling.
    """
    p_list = []
    # Generate a probability list along the kmer using profile matrix.
    for pos in range(len(sequence) - k + 1):
        p = 1
        kmer = sequence[pos: (pos + k)]
        # Calculate probability of each kmer by mutiplying all the p(letter)s.
        for col, letter in enumerate(kmer):
            p *= matrix[rows[letter]][col]
        p_list.append(p)
    index = rd.choices([i for i in range(len(sequence) - k + 1)],
                        weights=p_list, k=1)[0]  # 
    return sequence[index: index + k]

def ConsensusMotif(sequence_set):
    """Find the consensus motif for the input sequence set."""
    k = len(sequence_set[0])
    matrix = [[0] * k for i in range(4)]  # Initiate a zero-filled nested list.
    for col, colbases in enumerate(zip(*sequence_set)):
        for row, letter in enumerate('ACGT'):
            matrix[row][col] = colbases.count(letter)  # Count nuc in each col.
    convert = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    motif = ''.join([convert[col.index(max(col))] for col in zip(*matrix)])
    return motif
