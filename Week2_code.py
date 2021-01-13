#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 14:52:48 2020

@author: babemomo
"""
# %% exc. 1 Peculiar Statistics of the Forward and Reverse Half-Strands


def MinSkew(text):
    """
    Find a position in a genome where the skew diagram attains a minimum.

    Parameters
    ----------
    text : string
        A DNA string Genome.

    Returns
    -------
    All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).

    """
    text = str(text)
    skew = 0
    skewscore = {'A': 0, 'T': 0, 'C': -1, 'G': 1}
    walk = [0]
    minlist = []
    for i in text:
        skew = skew + skewscore[i]
        walk.append(skew)
    minskew = min(walk)
    for index, j in enumerate(walk):
        if j == minskew:
            minlist.append(index)
    return walk


# %% Exc.2 Hamming Distance Problem


def HammingDistance(p, q):
    """
    Compute the Hamming distance between two strings.

    We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi. For example, CGAAT and CGGAC have two mismatches. The number of mismatches between strings p and q is called the Hamming distance between these strings and is denoted HammingDistance(p, q).

    Parameters
    ----------
    Two strings of equal length.
    a : string
        A string of genome for comparison.
    b : string
        A string of genome for comparison.

    Returns
    -------
    The Hamming distance between these strings.

    """
    ham = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            ham += 1
    return ham


# %% Exc.3 Approximate Pattern Matching Problem


def ApproPatternMatch(pattern, text, d):
    """
    Find all approximate occurrences of a pattern in a string.

    We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.

    Strings Pattern and Text along with an integer d.

    Parameters
    ----------
    pattern : string
        k-mer for counting in the text.
    text : string
        A stretch of DNA sequence for pattern finding.
    d : integer
        Mismatch (Hamming distance) allowed for counting the pattern.

    Returns
    -------
    All starting positions where Pattern appears as a substring of Text with at most d mismatches.

    """
    l = len(pattern)
    mindex = []
    for i in range(len(text) - len(pattern) + 1):
        itext = text[i: (i + l)]
        if HammingDistance(pattern, itext) <= d:
            mindex.append(i)
    return mindex


# %% Exc.4 ApproximatePatternCount


def ApproPatternCount(pattern, text, d):
    """
    Strings Pattern and Text as well as an integer d.

    Parameters
    ----------
    pattern : string
        k-mer for counting in the text.
    text : string
        A stretch of DNA sequence for pattern finding.
    d : integer
        Mismatch (Hamming distance) allowed for counting the pattern.

    Returns
    -------
    count : Int
        The number of approximate match.

    """
    d = int(d)
    count = 0
    lpattern = len(pattern)
    for i in range(len(text) - lpattern + 1):
        itext = text[i: (i + lpattern)]
        if HammingDistance(pattern, itext) <= d:
            count = count + 1
    return count


# %% Exc.5 Frequent Words with Mismatches Problem
# temp = open('/Users/babemomo/Dropbox/My Mac (natekiMacBook-Pro.local)/Downloads/dataset_9_6.txt').read().split('\n')
# ApproPatternCount(temp[0], temp[1], temp[2])

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
            # try:
            #     neighbor_dict[j]
            # except KeyError:
            #     neighbor_dict[j] = ImmediateNeighbors(j)
            #     neighbor_set.update(neighbor_dict[j])
            #     temp_set.update(neighbor_dict[j])
    # # More intuitive code, however it has extra loop
    # neighbor_set = set([pattern])
    # neighbor_dict = {}
    # for i in range(d):
    #     count = iter(neighbor_set.copy())
    #     for j in count:
    #         neighbor_dict[j] = neighbor_dict.get(j, ImmediateNeighbors(j))
    #         neighbor_set.update(neighbor_dict[j])
    # return neighbor_set
    # with open("print_temp.txt", "w") as f:
    #     print(*neighbor_set, sep = '\n', file = f)
    # return neighbor_set


import time
def FrequentWordwithMismatch(text, k, d):
    start = time.process_time()
    thisdict = {}
    k_neighbor = set()
    for i in range(len(text) - k + 1):
        kmer = text[i: (i + k)]
        k_neighbor = k_neighbor.union(IterativeNeighbors(kmer, d))
    for kmer in k_neighbor:
        thisdict[kmer] = ApproPatternCount(kmer, text, d)
    maxcount = max(thisdict.values())
    maxlist = [i for i, j in thisdict.items() if j == maxcount]
    end = time.process_time()
    print("Time:", end - start)
    return maxlist


# %% Exc.6 

def ReverseComplement(Pattern):
    """Find the reverse complement of a DNA string."""
    thisdict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = [thisdict[i] for i in Pattern]
    # complement = []
    # for i in Pattern:
    #     complement.append(thisdict[i])
    return "".join(complement)[::-1]



def FrequentWordwithMismatch2(text, k, d):
    import time
    start = time.process_time()
    thisdict = {}
    k_neighbor = set()
    comp = ReverseComplement(text)
    for i in range(len(text) - k + 1):
        kmer = text[i: (i + k)]
        k_neighbor = k_neighbor.union(IterativeNeighbors(kmer, d))
    for i in range(len(comp) - k + 1):
        kmer = text[i: (i + k)]
        k_neighbor = k_neighbor.union(IterativeNeighbors(kmer, d))
    for kmer in k_neighbor:
        thisdict[kmer] = ApproPatternCount(kmer, text, d) + ApproPatternCount(kmer, comp, d)
    maxcount = max(thisdict.values())
    maxlist = [i for i in thisdict.keys() if thisdict[i] == maxcount]
    # for i in thisdict.keys():
    #     if thisdict[i] == maxcount:
    #         maxlist.append(i)
    end = time.process_time()
    print("Time:", end - start)
    return maxlist


if __name__ == '__main__':
    print(FrequentWordwithMismatch('ATGACCGGGATACTGATAGAAGAAAGGTTGGGGGCGTACACATTAGATAAACGTATGAAGTACGTTAGACTCGGCGCCGCCGACCCCTATTTTTTGAGCAGATTTAGTGACCTGGAAAAAAAATTTGAGTACAAAACTTTTCCGAATACAATAAAACGGCGGGATGAGTATCCCTGGGATGACTTAAAATAATGGAGTGGTGCTCTCCCGATTTTTGAATATGTAGGATCATTCGCCAGGGTCCGAGCTGAGAATTGGATGCAAAAAAAGGGATTGTCCACGCAATCGCGAACCAACGCGGACCCAAAGGCAAGACCGATAAAGGAGATCCCTTTTGCGGTAATGTGCCGGGAGGCTGGTTACGTAGGGAAGCCCTAACGGACTTAATATAATAAAGGAAGGGCTTATAGGTCAATCATGTTCTTGTGAATGGATTTAACAATAAGGGCTGGGACCGCTTGGCGCACCCAAATTCAGTGTGGGCGAGCGCAACGGTTTTGGCCCTTGTTAGAGGCCCCCGTATAAACAAGGAGGGCCAATTATGAGAGAGCTAATCTATCGCGTGCGTGTTCATAACTTGAGTTAAAAAATAGGGAGCCCTGGGGCACATACAAGAGGAGTCTTCCTTATCAGTTAATGCTGTATGACACTATGTATTGGCCCATTGGCTAAAAGCCCAACTTGACAAATGGAAGATAGAATCCTTGCATACTAAAAAGGAGCGGACCGAAAGGGAAGCTGGTGAGCAACGACAGATTCTTACGTGCATTAGCTCGCTTCCGGGGATCTAATAGCACGAAGCTTACTAAAAAGGAGCGGA', 15, 4))
# f = open('E_coli.txt', 'r')
# a = f.readline()
# FrequentWordwithMismatch2(a, 9, 1)
