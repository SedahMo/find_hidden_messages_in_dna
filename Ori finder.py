#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 13:38:03 2021

@author: Monica Ruan
"""
import time
start = time.process_time()
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    None

# %% Find the minimum G-C skew point.


def MinSkew(sequence):
    """Calculate the G-C skew scores along the sequence."""
    skew = 0
    skewscore = {'A': 0, 'T': 0, 'C': -1, 'G': 1}
    walk = [0]
    for i in sequence:
        skew = skew + skewscore[i]
        walk.append(skew)
    minskew = min(walk)
    minlist = [index for index, j in enumerate(walk) if j == minskew]
    return [walk, minlist]


# %% Find the Ori by check the coding window (+-250) nerby minimum skew point.


def HammingDistance(sequence_1, sequence_2):
    """Calculate hamming distance (sequence mismatch score)."""
    ham = 0
    for i in range(len(sequence_1)):
        if sequence_1[i] != sequence_2[i]:
            ham += 1
    return ham


def ApproPatternCount(pattern, sequence, d):
    """Calculate the approximate pattern in sequence with distance d."""
    d = int(d)
    count = 0
    lenth = len(pattern)
    for i in range(len(sequence) - lenth + 1):
        text = sequence[i: (i + lenth)]
        if HammingDistance(pattern, text) <= d:
            count = count + 1
    return count


def ImmediateNeighbors(pattern):
    """Generate a list of immediate neighbors of a pattern."""
    neighbors = [pattern]
    nuc = ['A', 'T', 'C', 'G']
    for i in range(len(pattern)):
        sym = pattern[i]
        for j in nuc:
            if sym != j:
                thislist = [pattern[:i], pattern[(i+1):]]
                neighbors.append(j.join(thislist))
    return neighbors


def IterativeNeighbors(pattern, d):
    """Generate a list of neighbors of a pattern with distance d."""
    neighbor_set = set([pattern])
    neighbor_dict = {}
    for i in range(d):
        count = iter(neighbor_set.copy())
        for j in count:
            neighbor_dict[j] = neighbor_dict.get(j, ImmediateNeighbors(j))
            neighbor_set.update(neighbor_dict[j])
    return neighbor_set


comdict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def ReverseComplement(Pattern):
    """Find the reverse complement of a DNA string."""
    complement = [comdict[i] for i in Pattern]
    return "".join(complement)[::-1]


def FreqWordwithApproMatch(sequence, k, d):
    """Find the approximate (with distance d) frequent k-mer in a sequence."""
    thisdict = {}
    k_neighbor = set()
    comp = ReverseComplement(sequence)
    # Gnerate all kmer neighbors in the sequence.
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i: (i + k)]
        k_neighbor = k_neighbor.union(IterativeNeighbors(kmer, d))
    # Generate all kmer neighbors in the complementary sequence.
    for i in range(len(comp) - k + 1):
        kmer = sequence[i: (i + k)]
        k_neighbor = k_neighbor.union(IterativeNeighbors(kmer, d))
    # Calculate the frequency of all the kmer neighbors and find the maximum.
    for kmer in k_neighbor:
        thisdict[kmer] = ApproPatternCount(kmer, sequence, d) + ApproPatternCount(kmer, comp, d)
    maxcount = max(thisdict.values())
    maxlist = [i for i in thisdict.keys() if thisdict[i] == maxcount]
    return maxlist


def main(file, output=True, header=False, k=9, d=1):
    """Define the main function."""
    # Open file and create a handle.
    with open(file, 'r') as handle:
        if header is True:
            next(handle)
            sequence = str()
            for i in handle:
                sequence = sequence + str(i.strip())
        else:
            sequence = handle.readline().strip()
    # Find minimum skewpoints and frequent words within the related window.
    skew_pt = MinSkew(sequence)
    worddict = {}
    for i in skew_pt[1]:
        window = sequence[(i - 250): (i + 250)]
        worddict[i] = FreqWordwithApproMatch(window, k, d)
    for m, n in worddict.items():
        print("Approximate frequent words at skew point",
              str(m) + ":",
              n)
    try:
        plt.plot(skew_pt[0])
        plt.show()
        if output is True:
            plt.plot(skew_pt[0])
            plt.savefig('Ori_output.png')
    except NameError:
        print("Matplotlib is not compatible with the interpreter.")


main('Salmonella_enterica.txt', header=True)
end = time.process_time()
print("Time:", end - start)
