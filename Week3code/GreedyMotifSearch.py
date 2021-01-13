#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 06:34:26 2021

@author: babemomo
"""
import numpy as np
import time
from Profile_most_probable_kmer import ProfileMostPKmer as findkmer

rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def ProfileMatrix(kmer_list, k, t, mode):  # nrow = t, ncol = k
    np_matrix = np.zeros((4, k))  # matrix.shape = 4, k
    for kmer in kmer_list:
        for col, letter in enumerate(kmer):
            np_matrix[rows[letter], col] += 1
    if mode == 'profile':
        return np_matrix / t
    if mode == 'score':
        return t * k - sum((np.amax(np_matrix, 0)))


def GreedyMotifSearch(dna_list, k, t):
    first_kmer_list = [sequence[0: k] for sequence in dna_list]
    best_matrix = first_kmer_list
    best_matrix_score = ProfileMatrix(first_kmer_list, k, t, mode='score')
    for pos in range(len(dna_list[0]) - k + 1):
        dna_list_iter = iter(dna_list)
        kmer_list = [next(dna_list_iter)[pos: pos + k]]  # Get first kmer.
        for sequence in dna_list_iter:
            promatrix = ProfileMatrix(kmer_list, k, len(kmer_list),
                                      mode='profile')
            prokmer = findkmer(sequence, k, promatrix)
            kmer_list.append(prokmer)
        kmer_list_score = ProfileMatrix(kmer_list, k, t, mode='score')
        # Compare the hamming scores, get the smaller one.
        if best_matrix_score > kmer_list_score:
            best_matrix = kmer_list
            best_matrix_score = kmer_list_score
    return best_matrix  # A np_matrix


if __name__ == "__main__":
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        params = line_list.pop(0).split()
        k = int(params[0])
        t = int(params[1])
        dna_list = line_list
    start = time.process_time()
    print(*GreedyMotifSearch(dna_list, k, t), sep='\n')
    end = time.process_time()
    print("Time:", end - start)
