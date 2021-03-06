#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 20:33:08 2021

@author: babemomo
"""
import random as rd
# import W4library as w4
import W4library_basic as w4  # without numpy
import time

# rd.seed(0)

def RandomizedMotifSearch(dna_list, k, t, m):
    """Generate consensus motifs by randomly choosing sequences from dna_list.

    Two loops. Outer one is for generating the initiative seeds, the result
    is super-sensitive to these initiatives. That's why in the inner iteration,
    we "jump out" immedeiately once the function is "derailed".
    """
    start = True
    for i in range(1000):
        ran_list = [rd.randint(0, len(dna_list[0]) - k) for i in range(t)]
        iter_seqs = [seq[ran_list[i]: ran_list[i] + k]
                     for i, seq in enumerate(dna_list)]
        if start:
            best_seqs = iter_seqs[:]
            start = False
        # Inner loop. Jump out of the loop if "derailed".
        while True:
            profile = w4.ProfileMatrix(iter_seqs, k, t, mode='profile')
            iter_seqs = [w4.ProfileMostPKmer(seq, k, profile)
                         for seq in dna_list]
            if w4.ProfileMatrix(iter_seqs, k, t, mode=m) < \
               w4.ProfileMatrix(best_seqs, k, t, mode=m):
                best_seqs = iter_seqs[:]
            else:
                break
    if m=='score':
        entroen.append(w4.ProfileMatrix(best_seqs, k, t, mode='score'))
        entroham.append(w4.ProfileMatrix(best_seqs, k, t, mode='ham'))
    if m=='ham':
        hamen.append(w4.ProfileMatrix(best_seqs, k, t, mode='score'))
        hamham.append(w4.ProfileMatrix(best_seqs, k, t, mode='ham'))
    print(w4.ProfileMatrix(best_seqs, k, t, mode='score'))
    return best_seqs


if __name__ == "__main__":
    with open('dataset_161_5.txt') as f:
    # with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        params = line_list.pop(0).split()
        k = int(params[0])
        t = int(params[1])
        dna_list = line_list
    start = time.process_time()
    # print(*RandomizedMotifSearch(dna_list, k, t), sep='\n')
    entroen=[]
    entroham=[]
    hamen=[]
    hamham=[]
    for i in range(100):
        rd.seed(i)
        print('entropy')
        RandomizedMotifSearch(dna_list, k, t, 'score')
        rd.seed(i)
        print('ham')
        RandomizedMotifSearch(dna_list, k, t, 'ham')
    end = time.process_time()
    print("Time:", end - start)
    print('entroen', sum(entroen)/len(entroen))
    print('hamen', sum(hamen)/len(hamen))
    print('entroham', sum(entroham)/len(entroham))
    print('hamham', sum(hamham)/len(hamham))
