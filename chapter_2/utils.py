import numpy as np
from typing import List

nucleotide_to_index_map = {
    "A": 0,
    "T": 1,
    "C": 2,
    "G": 3
}

index_to_nucleotide_map = "ATCG"

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

# TODO: make more efficient
def get_probability_of_k_mer(k_mer: str, profile_matrix: np.ndarray) -> float:
    prod = 1
    for (pos_index, nucleotide) in enumerate(k_mer):
        prod *= profile_matrix[pos_index, nucleotide_to_index_map[nucleotide]]
    return prod

    # This approach takes nearly 3 times as long as imperative solution!
    
    # probabilities = [profile_matrix[pos_index, nucleotide_to_index_map[nucleotide]] 
    #                  for (pos_index, nucleotide) in enumerate(k_mer)]
    # return np.prod(probabilities)

def iterate_kmers(dna: str, k: int) -> List[str]:
    return [dna[i:i+k] for i in range(len(dna) - k + 1)]

def iterate_kmers_tests():
    assert(iterate_kmers("ACTG", 2) == ["AC", "CT", "TG"])
    print("iterate_kmers tests passed")

def get_score(motifs: List[str], k: int) -> int:
    profile_matrix = get_profile_matrix(motifs, k)
    consensus_k_mer = consensus_from_profile_matrix(profile_matrix)
    np_consensus_k_mer = np.array(list(consensus_k_mer))
    return np.sum([np.sum(np.array(list(motif)) != np_consensus_k_mer) for motif in motifs])


def consensus_from_profile_matrix(profile_matrix: np.ndarray) -> str:
    return "".join([index_to_nucleotide_map[np.argmax(nuc_dist)] for nuc_dist in profile_matrix])

def run_tests():
    get_profile_matrix_tests()
    iterate_kmers_tests()
    print("== ALL UTILS TESTS PASSED ==")

run_tests()
