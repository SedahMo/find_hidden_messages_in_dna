#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 09:46:22 2021

@author: babemomo
"""
import random as rd
import W4library_basic as w4
import time


def GibbsSampler(dna_list, k, t, N, mode = 'v2', repeat=20):
    start = True
    for i in range(repeat):
        ran_list = [rd.randint(0, len(dna_list[0]) - k) for i in range(t)]
        iter_seqs = [seq[ran_list[i]: ran_list[i] + k]
                     for i, seq in enumerate(dna_list)]
        if start:
            best_seqs = iter_seqs[:]
            best_score = w4.ProfileMatrix(best_seqs, k, t, mode='ham')
            start = False
        # Inner loop.
        for j in range(N):
            index = rd.randint(0, t-1)
            # Generate the profile matrix with the rest of the temp_seqs.
            iter_seqs.pop(index)
            profile = w4.ProfileMatrix(iter_seqs, k, t-1, mode='profile')
            if mode == 'v1':
                iter_seqs.insert(index, w4.GibbsKmer(dna_list[index],
                                                     k, profile))
            if mode == 'v2':
                iter_seqs.insert(index, w4.ProfileMostPKmer(dna_list[index],
                                                            k, profile))
                # iter_seqs[index] = w4.ProfileMostPKmer(dna_list[index],
                #                                        k, profile)
            iter_score = w4.ProfileMatrix(iter_seqs, k, t, mode='ham')
            if iter_score < best_score:
                best_seqs = iter_seqs[:]
                best_score = iter_score
    print(w4.ProfileMatrix(best_seqs, k, t, mode='ham'))
    return best_seqs


## benchmark dataset
# a=[
# "ATGACCGGGATACTGATAGAAGAAAGGTTGGGGGCGTACACATTAGATAAACGTATGAAGTACGTTAGACTCGGCGCCGCCG",
# "ACCCCTATTTTTTGAGCAGATTTAGTGACCTGGAAAAAAAATTTGAGTACAAAACTTTTCCGAATACAATAAAACGGCGGGA",
# "TGAGTATCCCTGGGATGACTTAAAATAATGGAGTGGTGCTCTCCCGATTTTTGAATATGTAGGATCATTCGCCAGGGTCCGA",
# "GCTGAGAATTGGATGCAAAAAAAGGGATTGTCCACGCAATCGCGAACCAACGCGGACCCAAAGGCAAGACCGATAAAGGAGA",
# "TCCCTTTTGCGGTAATGTGCCGGGAGGCTGGTTACGTAGGGAAGCCCTAACGGACTTAATATAATAAAGGAAGGGCTTATAG",
# "GTCAATCATGTTCTTGTGAATGGATTTAACAATAAGGGCTGGGACCGCTTGGCGCACCCAAATTCAGTGTGGGCGAGCGCAA",
# "CGGTTTTGGCCCTTGTTAGAGGCCCCCGTATAAACAAGGAGGGCCAATTATGAGAGAGCTAATCTATCGCGTGCGTGTTCAT",
# "AACTTGAGTTAAAAAATAGGGAGCCCTGGGGCACATACAAGAGGAGTCTTCCTTATCAGTTAATGCTGTATGACACTATGTA",
# "TTGGCCCATTGGCTAAAAGCCCAACTTGACAAATGGAAGATAGAATCCTTGCATACTAAAAAGGAGCGGACCGAAAGGGAAG",
# "CTGGTGAGCAACGACAGATTCTTACGTGCATTAGCTCGCTTCCGGGGATCTAATAGCACGAAGCTTACTAAAAAGGAGCGGA",
# ]


# print(w4.ConsensusMotif(GibbsSampler(a, 15, 10, 2000)))


if __name__ == "__main__":
    with open('dataset_163_4.txt') as f:
    # with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        params = line_list.pop(0).split()
        k = int(params[0])
        t = int(params[1])
        N = int(params[2])
        dna_list = line_list
    start = time.process_time()
    # best_seqs = GibbsSampler(dna_list, k, t, 2000)
    seed_list = []
    seq_list = []
    score_list = []
    for i in rd.choices([i for i in range(100000000)], k = 30):
        # for N in range(2000, 2001):
            rd.seed(i)
            seed_list.append(i)
            best_seqs = GibbsSampler(dna_list, k, t, N)
            seq_list.append(best_seqs)
            score_list.append(w4.ProfileMatrix(best_seqs, k, t, mode='ham'))
    print('Best seed is ', seed_list[score_list.index(min(score_list))])
    print('Min score is ', min((score_list)))
    print('Average score is ', sum(score_list)/len(score_list))
    print('The scores are ', *score_list)
    print('Best sequence is ',
          *seq_list[score_list.index(min(score_list))], sep='\n')
    end = time.process_time()
    print("Time:", end - start)
