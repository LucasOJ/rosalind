from typing import List, Tuple
from math import inf
import numpy as np

from utils import get_profile_matrix, iterate_kmers, get_probability_of_k_mer, get_score

'''
dna_regions - list of t DNA nucleuodide strings
k - k-mer length
N - number of sampling iterations
'''
def gibbs_sampler(dna_regions: List[str], k: int, t: int, N: int) -> Tuple[List[str], int]:
    min_score = inf 

    assert(len(dna_regions) == t)
    max_k_mer_start_pos = len(dna_regions[0]) - k
    initial_motif_start_positions = np.random.randint(0, max_k_mer_start_pos + 1, t)
    motifs = [dna_regions[i][pos:pos+k] for (i, pos) in enumerate(initial_motif_start_positions)]
    
    best_motifs = motifs.copy()

    for _ in range(N):
        modify_dna_index = np.random.randint(t)
        # Make faster by copying then dropping one?
        motifs_excluding_selected = [motif for (i, motif) in enumerate(motifs) if i != modify_dna_index]
        profile_matrix = get_profile_matrix(motifs_excluding_selected, k)
        sampled_k_mer = profile_sample_k_mer(dna_regions[modify_dna_index], profile_matrix, k)
        motifs[modify_dna_index] = sampled_k_mer

        # Make one function
        score = get_score(motifs, k)

        if score < min_score:
            best_motifs = motifs.copy()
            min_score = score

    print("FINISHED SAMPLE")
    return (best_motifs, int(min_score))

def profile_sample_k_mer(dna_region: str, profile_matrix: np.ndarray, k: int) -> str:
    k_mers = iterate_kmers(dna_region, k)

    # i^th probability corresponds to i^th k-mer
    k_mer_probabilities = [get_probability_of_k_mer(k_mer, profile_matrix) for k_mer in k_mers]

    normalised_k_mer_probabilities = np.array(k_mer_probabilities) / sum(k_mer_probabilities)
 
    return np.random.choice(k_mers, p=normalised_k_mer_probabilities)

motifs = [
    "ACGAAATCAGCAATGCACCGACCAACGATATTCCCCCTCAAGCGGTACGGCTATCCCTTCATAGAACTATAAGAAGGCCTCAGAAACCCGTGCTCGGGTAAGACACGCGCCAAACGAGTTGGTGCGAAGCATGCATCAACGTAGCGCCGGGCTCCGGGCGCCGATGCTACCTACCTACCACGATCCCCCCTACCCGCAGTGAGGGGACAACGATTGAGAACACCTGTCGATCAAGTTCTCAATGTCAGTCTGGGGTAAACTGCTAATAAACGGGGGACACCTGTTATCGCATATCCCCTCGGACGAAATCAGCAATG",
    "CACCGACCAACGATATTCCCCCTCAAGCGGTACGGCTATCCCTTCATAGAACTATAAGAAGGCCTCAGAAACCCGTGCTCGGGTAAGACACGCGCCAAACGAGTTGGTGCGAAGCATGCATCAACGTAGTTGTCTTTAAACCCACGCCGGGCTCCGGGCGCCGATGCTACCTACCTACCACGATCCCCCCTACCCGCAGTGAGGGGACAACGATTGAGAACACCTGTCGATCAAGTTCTCAATGTCAGTCTGGGGTAAACTGCTAATAAACGGGGGACACCTGTTATCGCATATCCCCTCGGACGAAATCAGCAATG",
    "TTTTATGCGTTGACACTGTGTGGGAGGCTCATATACATCATGACATGACACCTAGGCCGGCTGACCGTGTCCGCACAGCACACTGTCAGCGGACTTTTTGTGGCTGGGACTGAGCTCAGCTTGAGTGACCTCTTGGTCGGCATCTGCAGAAGCCACCTAGGCAAACCAGTAGACCAGTGGGCCGGTACAGTGATAACATGATGTTGGGCGATAACCCAGGTCTCAAGACTTGCGTAGTCTCAAGTACAGTCAGTAACGCTACCCTAAAGCGGTACAACTGGTACAAGATAGCATTAATGGTTAGCGTTGCAATGTAG",
    "TGCGGAGCACTCCTATCGATCCTTGGGTCGCACCTGTATGTCCCGGTCCATCTGGCGATCGAAGCCGTAATCATCTGTCATTAAGCTCACTAGCGCATTGTTATCCGTCCGTTGAACTTCGAGGTCCAGCCGAGTGTACGCCATCGATTGGGGGAAAACCCAGTCGTGTTCGCCTAGACCGTGCCAGAGCCTGCCTGGAAAAGGTAAGTGAGACATTACTTGCATCTCTCGCAACCTCATGAATACCGTCCTTACAACAGGTCACAGAAGGCCCCGAAAACCCACTAAGCAATGTGACCTGGTTCCCGCGTCCCACG",
    "ATTTCCTAGTCGAGATTCTGCGACGCCGTCCAACGTTTTTGGTGTGCGTGGTCAGTAAGATCTTTGGCTTCAAAATCCCTGTAGTCGGTTATTAACGTTGGGCTTAAAGATAAATTGGCGGTGCCGAGCAGTGCGAGAGGGAAGGTCCCACGTGTTTGGAAGTCGTTCAATAGCGGCGCCGACTTTGATCATTTATACCAAACGATAGTAGGCGTTCAGACATGACGCACAAAGGCAGTACCGATAGAACCCCCATGCCCCAGGAACTTCACATACTCAGTTCAAATAGGAACGATGGCTCTGGAGGCTTTAGCGGG",
    "GCGTTGACCAGTGTCATTGGGCAGCGTGGCTATCCGGACGCATTCACTGGGTGATTGTTGGGGTCACTGCCACGCGTTAGTAATGTAAATAATTAATATGACACAACGTAGGCATGCGGGCTTAAACCCCTATCAAAAATGGAGAACGTGTCTATATTAGCGTATAATCGAGGGTATTCGAACCTCCATTCCAAGGGCCTTACAGTTTGGTACAGCCGGATGCAAGTATATGGCGTGGCTCCTACTGTAGCGACTTTGTTCCAGGAGTGATACACACTCCATTCCACTAGTGGAATGCAGACCTGAGCGTGCCAAGG",
    "GAGTAGCTATTGTCATAGCCACCCCCACGCCAGCGTCGGGCCCTGTCTACGGTCTTATAGTGTTCGCCAGTGGAAGATTTCGCGCTGGCAGGTCGCGAATTAGAGGGTGCTTGCCTTCTGTATTATGGAACGAGTAGGAAGTCGCCCGGTCGCCCGTATCACGCTATATAGTTTCAAAGTAATACCGGACACGAGTTAACATTAGCCCTGTGCGCTAATAACTTCGGCGCCAACCGGTCCTTGCCTGCCCCACCGTCATGGTGACGGTTTTGTGCTCTGGCTCTTCTGGGCTTAAACCACCTCTATGTGGACGCTCA",
    "AGTACAAGTCGTGGACATTGATTCTGGACAAAACAATCATATATATTATATACGGTAGTACCCGACATACTACGGCACGCATTACCTCAGAGATGGGACTTGGGCTTACTACCATGTGAGTTTCTCTCTGCTCAGCTGCGGCGAGTGAATTTATGGGGAAACAGCATTGAGTAAGCTATAGATACTATAAAAGTCGGCCACAACCATGTTCCCGTTGGCCCGTCTGGTCAGGGCCCGTCTGGTTCGCACGTTATTGGTCCTTGCACAGAGGATAGGGAGGAGGGACCTCATTTCGCATCACACACCCAAAGCCATAG",
    "CCAAGACTTGGGCTTAACTGCATACGAAGGGGACAATCGCAGCGGGCTTGCACGCCCATATGTTAGAGTATTACAGATCCAAAGCTCAAAAATGGACCCAGTATCTCAAACGAGCCGGGTGAGGGTAAGATCGGGTCGTGCCCGGGGCGTTTACGAGCTATAGATGCTATAAAGTAGTCCTGGTAATTGATGATGGGCAGTTGCACCTTCGGGTGGCCACAGCCTTCACTACCACGCTTTTCTACTTAAGACGAATGTAACTTGTAGGTCCCGGTCCCGGCTGTGACTCACTCAGTAACACGTCAAGCATTGCTCGA",
    "GCGACGGCCGCATTCGTCGTGGTAAGGACCTGCTTATAGGGCAGACGGTTGCTATGGGTATTCATCCCATGCTATCGGAAGACTGTCTAGGCATATATCGATGAACCTTAGCCTCCGATTGTTCAGGTACTTAGCACATTCGTGAGTAATAGTTGCCGCAAAGTCTTGGGGGGAAACCCACTGAGGTGCGGAAATTCGCTAGTCTACATCTGAACACGTAACGTTGACTCGGAATACGGCCCTGGTTGGGCGTTCGGCCAGGTAATATCGTGCGCCGAGGTAGCATGCTGTTTTTGCGATGCTCCTACCCCCATACC",
    "GCATTTCTCGCGAAATGCATCAATCTTGTTAGATTCTTGGCGCTAAACCCATTTGTGCTAAACAGACCACCCCTTTAAGGTGGGAGATCCATGGCCTTGAGCGCATCTACCCGTCTTCAACATAAGTGTCCGAGGGAACACATGCGTCAGATCCTCGTCGCGATCTTACTGCCATGGAATTTCACAGCTCGAATAACGAATCAAACCAGGACAATAAAATCATCAAGCGCGTTTATCCATGATCTCAATACAACACCAAAGCATCGCGCACGTCTCATCCCCGGCTTGCTCTTCTCAGTCTGTCCGCCACATTTACA",
    "AGAACTAGAGAGGCGCAACGACACGATCGAAGCCCTTAATTCACGTTCCAAAGTCTCCCTCATTCAATCTATTACCCTGGTAAAAATCGGATTAATCACTCAGTAGGAGCCTGGTGAGCCATGCGGGGTCAACGCTGGTTTTTACTTATTCTTAAACCCAGGAACTTAGCACATGAGCTCGAAACTTCGAGATGGATGGCTCTACATAGCCGCTATACTGGAGCGGTACCTTCTTATAGGGCGTGGGGGGTTTTAGGAGCATTGAGAACGTCTTCCAAGGACCAGATCGGCACCGCGAGGCCTACTAGAAGCACCGC",
    "TTGATAGACCCTGAAAAGGTAAATCAATCAGAGCCTACATGTAGCAATAGGCTGGTGAGATGACACACCGCATATGTGAAGGGTGGGCGCAACACGCGCCTTCCAAACGGTTGTGTGGCTTTTGGAAGCCCGAGGCAAGCCAACAGAAATATCCTGGCTACACTATCTCACCGATTCCAAACGCTAAGACTGAGCCCCTGAGGTCTTAATTTGCGCTATCCGGTTCAGCTGATAGTTGCTTTTAAACCCATTCTCCGCCATTGAAATACTTCTTATATAGTGCTATCCCAGTGGCCACCACACTCGTATAAGAGTTC",
    "CCTTTCGCTAGGCGCAGTTCGTTGTAACTAAGGTTAGTGGTCGATGTCCGAATGAATTTTGTATGCGAAACCCTTGCGTCATTAGTGTAGCTCGTGTTGGCCGGGACAGTAGTCACCCGCCCAAGCCCTGGTGAGTGGCCGGTGTACGGGCACCTGGGTACCATGTCGCTCCACCAGCCAGTTAGTCAAGCTTTGCATAGATGGGGGAGCAATTATCACCACCAACCTGGTGCCAAGGAAATGGATTTGGGCTGTGACCCACCATGTCGTTCAATCCCTATGGATGAGGCTATCACACACCGGGGATCCGACATCCT",
    "CCCGGTGTAAACCCTTCCGCATTGATGGGTTTCCTTGTTAGGACCCCGTGGTGCTTACAAAAGATGAAGTGCTGCAACTCTAAATTAGCTCTGTTACGTACTCGTCCGTCTATCCTGTTCACTACATTATTGCATGCGCCAACAAAGCCGAGTTTTCGATTATTTTACCCGGCTACATTGTAGCTAGTTGGGATATTACTATTTGGAGATAAACCCAACTAACTCACCAGCCTTCTCCGAACTCAAATACGAGAATGAGTCTAAGCTCCGGAGCTGTATGCGCATATGGTACTGCTGGCGATTGTGCATTCCCTCCA",
    "AACTGCGCTATTAGGCTGGTTCATGACACTTCACTCCAACTAATACTGTTAGAAGGTCCGGTTCGATATGCAGGACGCCAGTGGCGGCTTAAACCCATCGACTGACATTCGATTCCCATTGGGAGGAGGTAATAGTTAGTGGGAACGACCTTGTAGTAGCTGTTTACATCACTCGTGTGGTATGCGCTGATTACTCCACATTTTCCGCGACTAACATACGGCCAAATCATCCCGGGCTAGGGTGGCTCTCGGGTGGACCCCGACCATAACGAGACGGAGTCGCACCAAGCACATCGAGAATGGCAGACATCAAAAAC",
    "GCCACAGGGGGCTCTCTATCGCGTCTAGAGTAGCAATTCGATGGGCGCCATCGGGTCAAGTCGACTAGGTAGCACAGAATGTTTCCCCTCTGCTGCTTAAACCCACAGTAATGTAGCACTGGTTGTGATCGCAAAGGTTTGTATATGGTGATGAGACCCTGAAGAAGGTAGCTTAATCGAAGTCGAGCTAATTACATGTATTCAGGTTAGGTTCCGCACAGCCGCTTCGCTTGGGGGGGAATCTATTCACTTTGGATCTTCTGTCTAGCAGTATGCCGTGCACGAACTGGTCGTTGTGGATGTAACACGCAACCCAT",
    "TCGCTCTGCTTATAGAAGACGACTTAACTAGCCTTAGCCCCGGTGACATTCTAAACACAGCGGACTTTACCAGGGGGGAAGTCTGCCCAATAATCAAACAGTCGGCTTTATGTGATTACATACCCACTAGCGATCCACATCTTCATCATTTCAGTAGCAAAAGCTCCGGTTAGGGCTCGAGTACCCTCCACAATGCGGGGCCTCGAGGTGAGATCGTCCCGTACTTCGACCTTGGGCTTTGCCCCACTGTAACGCAGGCAAATGTAGCATCAGCTTATCACGGTCTTACCGCTTCTTCACTCTGAAGTATGCACCGA",
    "CAGGATTCATGTTTCAGTATGAGACCCCTCTGAGTTCTAGCCTGAATCTGAATGGCCAAACCTCGCAAGCTGGCGATGCCAATGGGACCCCGGCTCAAGATAATTCGAAAATACGAATCATTCCAATATACCAATGCCAATCCCGTGCCTTGAACCTTGCGCTTTCTTGGGCACCAACCCAAAAGACTCCCAATAGTCCCAAACCCGCACCAGGAAACCCTACCGATATTTGAGTCCATGGTCACGTAACTGAATTGGCTCAGTAGGCGCCTAGAGCATCCTCTTGGTTGTCGGCTGTGCACGAAGGCCTTGGGTCG",
    "CGACGACTATGCAAGAACTGGAATGCCATCCACGACTGCCTTCCGACGCATATGTCACGATTGCCATACTTCTAGTGTTAGTACCTTTCAAGTTTAGCATACGGCGATGAGGGATGTATGAACAGATCGACAATGGGTGGGTAAGTTGGCAAGGATACGAGATTTTCTCTCCTTCTGACACTTTGCGTCTGACACAGCCCCTTCCCTAAAGTTATTTTTGAGCTTTACCCGTTTGGGCTTAAACGACTCAGAACTCATGGGACACTTTCAACTCCATTGCCTTGCGTTACCTATATAATACCAGAGCACAGCCCTAA",
]

def trial():
    gibbs_sampler(motifs, 15, 20, 2000)

results = [gibbs_sampler(motifs, 15, 20, 1500) for _ in range(30)]
for k_mer in min(results, key=lambda x: x[1])[0]:
    print(k_mer)