#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 02:00:20 2021

@author: babemomo
"""

import time
import itertools as it
# from W3library import IterativeNeighbors as Neighbors
from DistanceBetweenPatternAndStrings import MinHamDistance as MinHam


def MedianString(k, dna_list):
    min_distance = None
    for pattern in [''.join(pattern) for pattern in it.product(
                                            ['A', 'C', 'T', 'G'], repeat=k)]:
        if min_distance is not None:
            d = MinHam(pattern, dna_list)
            if d < min_distance:
                motif = [pattern]
                min_distance = d
            elif d == min_distance:
                motif.append(pattern)
        else:
            min_distance = MinHam(pattern, dna_list)
            motif = [pattern]
    return motif


# def MediStrNeiboor(k, dna_list):
#     min_distance = None
#     for pattern in Neighbors('A' * k, k):
#         if min_distance is not None:
#             d = MinHam(pattern, dna_list)
#             if d < min_distance:
#                 motif = [pattern]
#                 min_distance = d
#             elif d == min_distance:
#                 motif.append(pattern)
#         else:
#             min_distance = MinHam(pattern, dna_list)
#             motif = [pattern]
#     return motif


if __name__ == "__main__":
    with open(input('Please input the path and press enter: \n')) as f:
        line_list = f.read().splitlines()
        k = int(line_list.pop(0))
        dna_list = line_list
    # start = time.process_time()
    # print('Neighbor generator', *MediStrNeiboor(k, dna_list))
    # end = time.process_time()
    # print("Time:", end - start)
    start = time.process_time()
    print('Itertool generator', *MedianString(k, dna_list))
    end = time.process_time()
    print("Time:", end - start)
