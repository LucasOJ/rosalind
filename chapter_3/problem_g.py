# https://rosalind.info/problems/ba3g/

from typing import List, Dict

def adj_list_to_adj_map(adj_list: List[str]) -> Dict[str, List[str]]:
    adj_map = {}
    for mapping_str in adj_list:
        [source, destinations_str] = mapping_str.split(" -> ")
        destinations = destinations_str.split(",")
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
def get_path_start_node(adj_map: Dict[str, List[str]]) -> str:
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
def traverse_from_node(adj_map: Dict[str, List[str]], start_node: str) -> List[str]:
    traversal = []
    stack = [start_node]
    while len(stack) > 0:
        # Get node at the top of the stack
        current_node = stack[-1]

        if current_node in adj_map:
            neighbourhood = adj_map[current_node]
            if len(neighbourhood) > 0:
                # If the current node has a neighbour left ...

                # ... remove the connecting edge from the graph ... 
                neighbour = neighbourhood.pop()
                
                # and visit the neighbour
                stack.append(neighbour)

                continue
        
        # If the current node has no neighbours left
        traversal.append(current_node)
        stack.pop()
    
    return traversal

def get_euler_path(adj_map: Dict[str, List[str]]):
    path_start_node = get_path_start_node(adj_map)
    return reversed(traverse_from_node(adj_map, path_start_node))

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
euler_path = get_euler_path(adj_map)
print("->".join(euler_path))


