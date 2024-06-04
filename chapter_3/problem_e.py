# https://rosalind.info/problems/ba3e/

from typing import List, Dict

def db_graph_from_k_mers(k_mers: List[str]) -> Dict[str, List[str]]:
    adj_map = {}
    for k_mer in k_mers:
        k_mer_prefix = k_mer[:-1]
        k_mer_suffix = k_mer[1:]
        if k_mer_prefix in adj_map:
            adj_map[k_mer_prefix].append(k_mer_suffix)
        else:
            adj_map[k_mer_prefix] = [k_mer_suffix]
    return adj_map

k_mers = [
    "GAGG",
    "CAGG",
    "GGGG",
    "GGGA",
    "CAGG",
    "AGGG",
    "GGAG",
]

db_graph = db_graph_from_k_mers(k_mers)

for (k, v) in db_graph.items():
    print(f'{k} -> {",".join(v)}')
