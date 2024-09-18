#!/usr/bin/env python3
import sys

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']
edges = {}
k = 3

# Dictionaries use sets for keys, and keys can not be duplicates, so we use dictionaries to create the kmer unique sequence libraries
for read in reads:
    for i in range(len(read) - k):
        kmer1 = read[i: i+k]
        kmer2 = read[i+1: i+1+k]
        edges.setdefault((kmer1, kmer2), 0)
        # edges[(kmer1, kmer2)] += 1 --> Is not required

# Used to save the digraph and apply it in the python code for graphviz
print('digraph {')
for kmer1, kmer2 in edges:
    print(f"{kmer1} -> {kmer2}")
print('}')