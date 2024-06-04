# https://rosalind.info/problems/ba3g/

from typing import List, Dict, Set

def adj_list_to_adj_map(adj_list: List[str]) -> Dict[str, Set[str]]:
    adj_map = {}
    for mapping_str in adj_list:
        [source, destinations_str] = mapping_str.split(" -> ")
        destinations = set(destinations_str.split(","))
        adj_map[source] = destinations
    return adj_map

# Adds default degree statistics if none have been recorded yet for the given node
def try_add_default_stats_for_node(node: str, degree_stats: Dict[str, Dict]):
    if not (node in degree_stats):
        degree_stats[node] = {
            "in_degree": 0,
            "out_degree": 0
        }

# Looks for unbalanced nodes
# Assumes unique unbalanced "start" and "end" nodes
def get_path_start_node(adj_map: Dict[str, Set[str]]) -> str:
    degree_stats_map = {}
    for (source, destinations) in adj_map.items():
        try_add_default_stats_for_node(source, degree_stats_map)
        degree_stats_map[source]["out_degree"] += len(destinations)

        for destination in destinations:
            try_add_default_stats_for_node(destination, degree_stats_map)
            degree_stats_map[destination]["in_degree"] += 1
    
    for (node, degree_stats) in degree_stats_map.items():
        if degree_stats["out_degree"] > degree_stats["in_degree"]:
            return node
    raise Exception("No nodes with greater out degree than in degree")

# Returns reversed Euler Path starting from a given node in the graph specified by the given adjacency_map
# Once a node has been traversed, it will have no neighbours left
def traverse_from_node(adj_map: Dict[str, Set[str]], node: str) -> List[str]:
    # In case the node has no neighbours
    if not (node in adj_map):
        return [node]
    
    traversal = []

    neighbourhood = adj_map[node]

    # Traverse all neighbours until there are none left
    while len(neighbourhood) > 0:
        neighbour = neighbourhood.pop()
        # This traversal may remove other direct neighbours of `node`
        adj_traversal = traverse_from_node(adj_map, neighbour)
        traversal.extend(adj_traversal)
        
    traversal.append(node)
    return traversal


adj_list = [
    "0 -> 2",
    "1 -> 3",
    "2 -> 1",
    "3 -> 0,4",
    "6 -> 3,7",
    "7 -> 8",
    "8 -> 9",
    "9 -> 6",
]

adj_map = adj_list_to_adj_map(adj_list)
path_start_node = get_path_start_node(adj_map)
traversal = reversed(traverse_from_node(adj_map, path_start_node))
print("->".join(traversal))


