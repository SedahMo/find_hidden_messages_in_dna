#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 13:38:03 2021

@author: Monica Ruan
"""
import time
import argparse
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


# %% Find the DnaA box by check the OriC window nerby minimum skew point.


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
    """Generate the reverse complement of a DNA string."""
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


def main(file, mode='center', output=True, header=False, k=9, d=1, length=1000,
         intersection=False):
    """Find minimum C-G skew points and frequent words.

    Works for Escherichia coli and Salmonella enterica.
    The default length of the window around skew point is length=500,
    kmer length kmer=9, mismatch tolerance d=1.
    Mode can be set to "start", "center" or "back" which is the position of the
    window related to the first minimum skew position.
    Mathcuinon controls whether using the union of the results from d=1 and d=2.
    """
    # Open file, create a handle, parse the text and make the sequence.
    with open(file, 'r') as handle:
        if header is True:
            line_list = handle.read().splitlines()[1:]
            sequence = ''.join(line_list)
        elif header is False:
            line_list = handle.read().splitlines()
            sequence = ''.join(line_list)
    # Print the specification of the call
    print("The mode is:", mode,
          "\nThe window length is:", length,
          "\nThe k-mer length is:", k)
    if intersection is True:
        print("The mismatch tolerance is equal or less than", str(d) + ";",
              "results are shown as the intersection of all tolerance under",
              d)
    if intersection is False:
        print("The mismatch tolerance is:", d)
    # Find minimum skewpoints and frequent words within the related window.
    skew_pts = MinSkew(sequence)
    print("The total length of the genome is", len(sequence),
          "\nMinimum skew points are located at:", *skew_pts[1])
    worddict = {}
    first_skew_pt = skew_pts[1][0]
    if mode == 'start':
        window = sequence[(first_skew_pt): (first_skew_pt + length)]
    elif mode == 'center':
        window = sequence[(first_skew_pt - int(length/2)): (first_skew_pt + int(length/2))]
    elif mode == 'end':
        window = sequence[(first_skew_pt - length): (first_skew_pt)]
    if intersection is False:
        worddict[first_skew_pt] = FreqWordwithApproMatch(window, k, d)
    elif intersection is True:
        set_words = set(FreqWordwithApproMatch(window, k, d=0))
        if d >= 1:
            for i in range(1, d+1):
                set_words &= (set(FreqWordwithApproMatch(window, k, d=i)))
        worddict[first_skew_pt] = list(set_words)
    # Check if the specified DnaA box is found.
    DnaA_box = set(["TTATCCACA", "TTATCCAAA",
                   "TTATCAACA", "TTATCAAAA"])
    if set(list(worddict.values())[0]) & DnaA_box != set():
        answer = "Yes"
    else:
        answer = "No"
    for m, n in worddict.items():
        print("Find total " + str(len(list(worddict.values())[0])),
              "approximate frequent words within",
              str(length) + 'bps',
              str(mode) + "ed", "at",
              str(m),
              "(..." + str(sequence[m: (m+10)]) + "...)",
              ":\n",
              n,
              "\nFind the defined DnaA box (5'-TTATC[CA]A[CA]A-3')?",
              *[answer])
    try:
        end = time.process_time()
        print("Run time:", end - start)
        plt.plot(skew_pts[0])
        plt.show()
        if output is True:
            plt.plot(skew_pts[0])
            plt.savefig('Ori_output.png')
    except NameError:
        print("Matplotlib is not compatible with the interpreter.")


try:
    if __name__ == '__main__':
        # Add shell command parser.
        parser = argparse.ArgumentParser()
        parser.add_argument('--file')
        parser.add_argument('--mode')
        parser.add_argument('--output')
        parser.add_argument('--header')
        parser.add_argument('--k')
        parser.add_argument('--d')
        parser.add_argument('--length')
        parser.add_argument('--intersection')
        args = parser.parse_args()
        main(file=str(args.file) if args.file is not None else 'Salmonella_enterica.txt',
             mode=str(args.mode) if args.mode is not None else 'center',
             # Argument following else is the default.
             output=False if args.output in ['False', 'F'] else True,
             header=False if args.header in ['False', 'F'] else True,
             k=int(args.k) if args.k is not None else 9,
             d=int(args.d) if args.d is not None else 1,
             length=int(args.length) if args.length is not None else 1000,
             intersection=True if args.intersection in ['True', 'T'] else False)
except:
    print("===========================Error occured. Did you " +
          "correctly specify header?===========================\n",
          "===========================Use --head True for files with header" +
          "===========================")
    raise

# f = ''.join(open(file, 'r').readlines()[1:])
# TTATCCACA
# https://string-db.org/network/220341.16505650
# https://string-db.org/network/316407.85676342
# (dnaA box): 5'-TTATC[CA]A[CA]A-3'

