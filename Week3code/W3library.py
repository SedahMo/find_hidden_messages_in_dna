#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 15:27:01 2021

@author: babemomo
"""

def HammingDistance(p, q):
    ham = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            ham += 1
    return ham


def ApproPatternCount(pattern, text, d):
    d = int(d)
    count = 0
    lpattern = len(pattern)
    for i in range(len(text) - lpattern + 1):
        itext = text[i: (i + lpattern)]
        if HammingDistance(pattern, itext) <= d:
            count = count + 1
    return count


def ImmediateNeighbors(pattern):
    neighborhood = [pattern]
    nuc = ['A', 'T', 'C', 'G']
    for i in range(len(pattern)):
        sym = pattern[i]
        for j in nuc:
            if sym != j:
                thislist = [pattern[:i], pattern[(i+1):]]
                neighborhood.append(j.join(thislist))
    return neighborhood


def IterativeNeighbors(pattern, d):
    neighbor_set = set([pattern])
    temp_set = set([pattern])
    neighbor_dict = {}
    for i in range(d):
        if i == 0:
            count = iter(set([pattern]))
        else:
            count = iter(temp_set - (neighbor_set & temp_set))
            neighbor_set |= temp_set
            temp_set = set()
        for j in count:
            neighbor_dict[j] = ImmediateNeighbors(j)
            temp_set.update(neighbor_dict[j])
        if i == d - 1:
            neighbor_set |= temp_set
    return neighbor_set


def FrequentWordwithMismatch(text, k, d):
    thisdict = {}
    k_neighbor = set()
    for i in range(len(text) - k + 1):
        kmer = text[i: (i + k)]
        k_neighbor |= IterativeNeighbors(kmer, d)
    for kmer in k_neighbor:
        thisdict[kmer] = ApproPatternCount(kmer, text, d)
    maxcount = max(thisdict.values())
    maxlist = [i for i, j in thisdict.items() if j == maxcount]
    return maxlist