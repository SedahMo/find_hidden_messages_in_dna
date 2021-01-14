#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 18:30:27 2021

@author: babemomo
"""
import time

def greedymotifsearch(dna, k, t):
    bestmotifs=[string[:k]for string in dna]
    for j in range(len(dna[0])-k+1):
        motifs=[dna[0][j:j+k]]
        for i in range(1,t):
            motifs+=[profilemostprobablekmer(dna[i],k,profile(motifs))]
        if score(motifs) < score(bestmotifs):
            bestmotifs=motifs
    return bestmotifs

def profile(motifs):
    transposed=[list(row) for row in zip(*motifs)]
    n=len(motifs)
    profile={nucleotide:[i.count(nucleotide)/n for i in transposed] for nucleotide in 'ACGT'}
    return profile   

def score(motifs):
    transposed=[list(row) for row in zip(*motifs)]
    counted=[[i.count(nucleotide) for nucleotide in 'ACGT'] for i in transposed]
    scored=sum([len(motifs)-max(i) for i in counted])
    return scored

def profilemostprobablekmer(text,k,profile):
    probability=[]
    for i in range(len(text)-k+1):
        compute=1
        for j in range(k):
            compute=compute*(profile[text[i+j]][j])
        probability.append(compute)
    idx=probability.index(max(probability))
    return text[idx:idx+k]


if __name__ == "__main__":
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        params = line_list.pop(0).split()
        k = int(params[0])
        t = int(params[1])
        dna = line_list
    start = time.process_time()
    print(*greedymotifsearch(dna, k, t), sep='\n')
    end = time.process_time()
    print("Time:", end - start)
    
 