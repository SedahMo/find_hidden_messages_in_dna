#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 04:53:23 2021

@author: babemomo
"""
import numpy as np
import time

rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def ProfileMostPKmer(sequence, k, np_matrix):
    max_p = 0
    motif = []
    for i in range(len(sequence) - k + 1):
        p = 1
        kmer = sequence[i: (i + k)]
        for col in range(k):
            letter = kmer[col]
            p *= np_matrix[rows[letter], col]
        if p > max_p:
            max_p = p
            motif = [kmer]
        elif p == max_p:
            motif.append(kmer)
    return motif


if __name__ == "__main__":
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        sequence = line_list.pop(0)
        k = int(line_list.pop(0))
        pro_matrix = [list(map(float, row.split(' '))) for row in line_list]
        np_matrix = np.array(pro_matrix)
    start = time.process_time()
    print(*ProfileMostPKmer(sequence, k, np_matrix))
    end = time.process_time()
    print("Time:", end - start)
