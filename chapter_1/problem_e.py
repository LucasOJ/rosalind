# https://rosalind.info/problems/ba1e/

from typing import List

# Finding all k-mer positions, and indentifying which correspond to (L, t)-clumps

def find_clumped_kmers(genome: str, k: int, L: int, t: int) -> List[str]:
    kmer_map = {}
    for i in range(len(genome) - k):
        kmer = genome[i:i+k]
        if not(kmer in kmer_map):
            kmer_map[kmer] = []
        kmer_map[kmer].append(i)

    # kmer_map is large: The keys alone are O(|genome|k)
    # Could make smaller by storing just the index of the first location

    # Parallelisable with map-reduce
    clumped_kmers = [kmer for (kmer, locs) in kmer_map.items() if is_clump(locs, k, L, t)]

    return clumped_kmers
            
# Gets more expensive the more k-mer repetitions we see
def is_clump(locs: List[int], k: int, L: int, t: int) -> bool:

    # Sliding window of t adjacent k-merk positions
    for i in range(len(locs) - t + 1):
        # Add k because locs[i] is the start position of the i th k-mer
        t_interval_end = locs[i + (t-1)] + k 
        t_interval_start = locs[i]
        t_interval_width = t_interval_end - t_interval_start

        # If true, then t adjacent k-mers in some window of length L
        if t_interval_width <= L:
            return True

    return False

# Keeping track of k-mer frequencies in a sliding L-window

# TODO