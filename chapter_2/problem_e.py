import numpy as np
from typing import List

nucleotide_to_index_map = {
    "A": 0,
    "T": 1,
    "C": 2,
    "G": 3
}

index_to_nucleotide_map = "ATCG"

def greedy_motif_search(dna_regions: List[str], k: int, t: int) -> str:
    assert(len(dna_regions) == t)
    
    consensus_pairs = []

    # For all k-mers in the first dna_region...
    for inital_motif in iterate_kmers(dna_regions[0], k):
        motifs = [inital_motif]

        # ...one by one extend with the most likely k-mer from each of the dna_regions
        for dna_region in dna_regions[1:]:
            profile_matrix = get_profile_matrix(motifs, k)
            most_likely_k_mer = get_most_likely_k_mer(dna_region, profile_matrix, k)
            motifs.append(most_likely_k_mer)

        profile_matrix = get_profile_matrix(motifs, k)
        consensus_k_mer = consensus_from_profile_matrix(profile_matrix)
        score = get_score(motifs, consensus_k_mer)
        consensus_pairs.append((score, motifs, consensus_k_mer, profile_matrix))
    
    # Choose motifs with lowest score
    min_score_pair = min(consensus_pairs, key=lambda x: x[0])
    return min_score_pair[1]

# Takes a list of k-mers each of length k
# Returns the profile matrix of the k-mers
def get_profile_matrix(k_mers: List[str], k:int) -> np.ndarray:
    for k_mer in k_mers:
        assert(len(k_mer) == k)

    nucleotide_counts = np.zeros((k, 4)) # 4 is alphabet size {A,T,C,G}

    for k_mer in k_mers:
        for (nucleotide_pos, nucleotide) in enumerate(k_mer):
            nucleotide_index = nucleotide_to_index_map[nucleotide]
            nucleotide_counts[nucleotide_pos, nucleotide_index] += 1

    # Laplace's rule of succession to get rid of 0s
    nucleotide_counts += 1

    # Each k-mer contributes 1 count at each position, 
    # and we added 4 extra count in Laplace's rule of succession
    nucleotide_count_sum = len(k_mers) + 4

    profile_matrix = nucleotide_counts / nucleotide_count_sum

    return profile_matrix

def consensus_from_profile_matrix(profile_matrix: np.ndarray) -> str:
    return "".join([index_to_nucleotide_map[np.argmax(nuc_dist)] for nuc_dist in profile_matrix])

def get_most_likely_k_mer(dna_region: str, profile_matrix: np.ndarray, k: int) -> str: 
    k_mer_probabilities = [get_probability_of_k_mer(k_mer, profile_matrix) 
                           for k_mer in iterate_kmers(dna_region, k)]
    max_prob_k_mer_index = np.argmax(k_mer_probabilities)
    return dna_region[max_prob_k_mer_index:max_prob_k_mer_index + k]

def get_probability_of_k_mer(k_mer: str, profile_matrix: np.ndarray) -> float:
    probabilities = [profile_matrix[pos_index, nucleotide_to_index_map[nucleotide]] 
                     for (pos_index, nucleotide) in enumerate(k_mer)]
    return np.prod(probabilities)

def get_score(motifs: List[str], consensus_k_mer: str) -> int:
    np_consensus_k_mer = np.array(list(consensus_k_mer))
    return np.sum([np.sum(np.array(list(motif)) != np_consensus_k_mer) for motif in motifs])

def iterate_kmers(dna: str, k: int) -> List[str]:
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]

### TESTS ###

def get_score_tests():
    assert(get_score(["abc"],"abc") == 0)
    assert(get_score(["abc"],"aaa") == 2)
    assert(get_score(["abc","daa"], "aaa") == 3)
    print("get_score tests passed")

def get_profile_matrix_tests():
    challenge_1 = get_profile_matrix(["ACGTA"], 5)
    expected_1 = np.array([
        [2, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 1, 1, 2],
        [1, 2, 1, 1],
        [2, 1, 1, 1],
    ]) / 5
    assert(np.array_equal(challenge_1, expected_1))

    challenge_2 = get_profile_matrix(["AAA","TTT"], 3)
    expected_2 = np.array([
        [2, 2, 1, 1],
        [2, 2, 1, 1],
        [2, 2, 1, 1],
    ]) / 6
    assert(np.array_equal(challenge_2, expected_2))

    print("get_profile_matrix tests passed")

def iterate_kmers_tests():
    assert(iterate_kmers("ACTG", 2) == ["AC", "CT", "TG"])
    print("iterate_kmers tests passed")

def greedy_motify_search_tests():
    challenge = greedy_motif_search([
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG",
    ], 3, 5)
    expected = [
        "TTC",
        "ATC",
        "TTC",
        "ATC",
        "TTC",
    ]
    assert(np.array_equal(challenge, expected))
    print("greedy_motify_search tests passed")

def run_tests():
    get_score_tests()
    get_profile_matrix_tests()
    iterate_kmers_tests()
    greedy_motify_search_tests()
    print("== ALL TESTS PASSED ==")

# run_tests()

### ROSALIND PROBLEM ###

motifs = [
    "AGTCATCCGCTTCTGTATATGCACATATGTCATAAAGAATTTGTTGAGCTCCGGATTCGTCGTTTCTGTTAATATCGCGCACAGGTTCAAGTTCTGCGCGAGAGACTCATATAGACACAAGCGTAGTGTAGGGCTCATATGAGGTCACCAAATTGC",
    "CCCCAATTCTGGTCCACCGTCACTGGCCCGATAAGCCGACGCGCCGGCCTCCGGATAAGTTGCGGTAAACTACTCCTCGGATCGCGTCGGAGTCCGGCTTGCGGACATCCCGAAATGTTTTGTAGCCGTCCGAGGTGGCCGCCTGTATTAACTTCC",
    "TTCGCCTGAATTACGCGCTTAGGATACTAAGAGAGCTTATACGCGTGTAACCTTTAAGCTCTCCCGATTGGGCATTGGACAGATTTGAACCCTTGCAGTGATGCTGTCGCGGGTCGGAACGAAACCTCGTCGCTCGGTGTAGAGTAGACGTGATAG",
    "CCCCCTCGGCGGCTCCCGATGAGCTTCCACCAATTCCTGTAAATGCGAGTGCCCTAGACTCTAGATCTGCAACGGCACGATTTGAAACGGCGATTGCGAGAAGGGCAAACTGCAACTGCTTATAACGAATCACGGGGAAAAATAGCTGAACTATGC",
    "GATGGAGTCCACCAGTGCGGGGAAGCGTAAGCGCTACCCAAAGCGCCAGCTATGAAATAAGACTAGAGGATATCAAAACCGGACACTTTGCTCATGCTCCAGATCAGAGTTATTCTTGCTCTAGCATTGGACCGAGTCCCACCCTTTGGATCTCAT",
    "CTGCAGTGGGAACTCCAGATCGGCACTTGTGCCAGATGATCGCATGCCTGCTACGTATACGTATCGGGAGTATACGGGAGTCTGTTATGCGCTAAAAAGAATTAGATCGATCCCGAACCTATTACTACAGGGGAGCAAACTACGACGAGCAGTGGG",
    "CCTTCTGCGGTACCTTACATGCTCTATACGGGAAAATGTTAAGTACCTGCAATCATTTCAAATCCTCCCCCTTGTCCGTCGCATTAACTCAGGGCCATGCCCCTGATCATCCGGTCTACACTCCTGATCTGCGTTACTATAACATAGCACTCTACG",
    "ATACGTGAATCACGCGATGCTCTCTATGGGACAAGTAGGCGCAATACAGATCAGTACTTACTCCAGATGCGGTTGAATTCAACTCGCTAGTGCTAGTTTCGGAAAATAGTCGTACAGGAACTGACGTAATTGAGTGAACCCTTTGCAGAGATGTTA",
    "ATGGGTAGTAGCATTAACTTCTGCGATACAAGATACGCATCCCGCAAGAGCAGCAGTTACGGATGTGGCCGTATTCTGTGAGGATTATCTCGTGTCGATGGCAAGTGTAAGTACGGTATTTCGGGGCGCCCGCTCCGGATTAGCGCATAGACTTTA",
    "GGGCAGCAAGCTATGCCCCTTCTTTGGATGTGTAAGCTCCCGATCTGCTCGTAGTTTAAGTCGGCCTCGTTTTTTTACGCCACTCGAGGGTCAATTCAGTACGATTGTTTCATTAGCATTAACTCTGCAAGAAGCGCACCGAATTACATAAGTGCT",
    "CCGGTTACGCATCGTCTCCCAGTTTCCACGGCGTTTCGCAAAAATCCGATCGTTGCAAGCGTTAAGATATGAGCGTGTAGTTTGTGCCAACTGGGAGTTATCTGTCCTCTCGAGGGCGACGACTCTGGACTAAATACCCCGTCTCTCCTGATTGGA",
    "AGTCTCACGAATTCTACTAGGGCGTGTCATCCTCCTCAGAACAGCCGCCATTATCCTGACGTTTAATTTATAGACTCATTCGTTCCGAGCTAATACCGAGATCATCGCTGGTGTTACCCGCTCCCGATTTGCCCACATATGGCTAGTCCGTCGAGG",
    "CTCCTGATGGGTGACGCCTGTTTTCCGAGTAAGACCCAAAACCGCCACGCGATCAGTCTTGGCCTTGGGTCAGTACGGAGTTGTGTTGTATCCCGTCTCAGGACTGCAAGCTTTCTTCTTGTCTTTGCAGTTCTTATGTGGATACGTTCCCCTTAA",
    "TGACCTAAACTGTATAAGCCAGGCCGGTGATCACTGGGGAACTCGTTTGAGCCAGCAGCCATGTATGTGCCCACCTACCTGCGAAGCTGGTGGGCTATGGTTGGGCAGCTCCCGATCCGGGACTTTTTATGAGCGATTACTATGGATTTAGACTGT",
    "CCTCTGGCATATTAATGCTTCGGACGCGGCGGGACCTCGGTTTCTTATAGAAATGGCGGCCACGTCCTGCCATTCGGTCATTCACGTGAGATTGTATCATTGACGGATTCTTCAATTAGTCTCCAGATGCGATGCATCTACTTGGCGGGGAACGGC",
    "GTCCTTCTTTAATCTATTTACCTATAAGCTCGAGTTAAGAAGCTTTGTATCCATTGGCATCTCCGGATAAGAAGTTTTACAAGTAACTACGTCATGGCCTGAGCGGGAAAGGTGGTCATGTGCAGTCAATTCGAGGAAGGTCCGTCATCCGCTCGG",
    "AACGATGTATGCCTCCTGATAGGCCGCTCGTAACATAACTCCCTATTTGGGACTGTAGCCGGATCAAGGTTGCTTATTTCTATCCACTGCACTTCATTCTGTTCCCCAGGTTAGTTCGCTAATATTGGGTGGCACGATTGATATCTCCTGAACATT",
    "CGGGATGCTGGCGCCACGAGTGCTTACACATTTGATATATGACACAGGTATTTCGCTACCCGTTGCACTAGCCCTTCTTTCACAATTTACGTTCCTGACACGACTAGTTTTGAAGATCTTCTCCGGATGTGACGCAACTGCGGCTACGTTACATCA",
    "TCACTTGAAGCGCCCCGGTACTGACACGAGGAACAGTGAAACATACAACTCCTGATTCGGTGAGCCTATGCCAATACGAGTGCCACTACCACGTTAGACGTGACTTAGCGTGGTACCGCGCTCGCGATAGTCTAGCGGTCAAGCCGTTGAGCTATG",
    "ACGCTTTAGCCTTGAGGGACTACGACGAATAACTACTAAATCACATAGCTCCAGATCTGGTCGAATGCCAGCGAGGACCTGGACAAACGGCTTCGAGGCCCCTACTGTGCGGGATCATTAAGACGTTAAAACTCACTGAAAATATTTGGCCACCTC",
    "CAGCAGGGTTTTTAAGTGTCACATCAGATCACTGTAACAGTGCCAAAGTTACTTTTTCCGCTCCAGATGCGTCTAACGGCCGCATGTCCAGTGATAGCAGGGCCGCCTTGGCATGGGCCATCTAGGACTCGCCTGGACGTTGCCTCCGACTGGCTG",
    "CTCCTGATGGGATCTGGCATTTAGAAATGTACCTGTAAACGGAAAAACTCTCGTATGGGGGTTCGCACTTTTTCTATGCGTGGCCCCGAGAGTGAGGGTTAAACCACGGAGGTGTTCCCAGAGATGTCATCCGTGTCTTTGCAGTGAGTAAGGCTA",
    "GCGGTCGTCCTTAGCGTTCCCACATATACAGGCGTTATAGGCCCAGGTACCAGGCTTTTGTCCTCGTGAGTGCATGAGCCTTGTCTCCCGATTCGGGTAGCAAAGCGATGAAAAACCCCTGAGGTGGTGGAGTGCTTACGTTTATTCATTCGTAAG",
    "TGGGGACCCTAAGATAGCTGGCTCTGTTAATTACCAATCCATTAAATCCGCATCTGTCCGCTCCGCTCGTTAGTTAGTAGTGTTCGACTTTATCACCTCCTGATTTGCCACAAGGTGGGCATTCGATACCGGCGCAAGGTGGCAGCTTTTACCTCA",
    "TTGCCATTGTTGACTTCATTATCGGTGCTCCTCAAATCGGAGGGGGCCCTCCAGATTAGAAAGCCGTGGCTACCCGATATTTCCGCTCGCAGCCCGCTTGCCCCTTGTGAAGCTGGCGACGAGCCGCCGATAACTACCCTTGGCCCCTTCCGAATA",
]

for motif in greedy_motif_search(motifs, 12, 25):
    print(motif)
