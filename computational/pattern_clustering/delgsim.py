#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 10:41:10 2025

@author: lucaci
"""

"""
This code contains the functions needed to compute the graph similarity measure desciber in Champin and Solnon (2003).
Given two graphs G_i = (V_i, E_i), i = 1,2, a relation m subset V_1 x V_2:
    
    sim(G_1, G_2) = max_{m} score(m, G_1, G_2)/f(G_1 intersect(m) G_2)
    
    
    score(m, G_1, G_2) = f(G_1 intersect(m) G_2) - g(splits(m))
    
    G_1 intersect(m) G_2 = {vertices in G_1 or G_2 which are related to vertices with the same label}
                                U
                            {edges in G_1 or G_2 whose endpoints are related to vertices connected by an edge with the same label}
                            
    splits(m) = {(v, m(v)): v in V_1 U V_2, |m(v)| > 1}
    
    f(F) = sum_{v in F}w(v) + sum_{e in F}w(e)
    
    g(splits(m)) = sum_{v in V_1 U V_2: |m(v)| > 1} w(v, m(v))
    

In this script the following objects are considered:
    - a relation m is described as a dictionary,
        in which each entry is a vertex/edge of V_1 U E_1 U V_2 U E_2 associated with a list of nodes/edges
        that are related to it via m. It is assumed that every update of m preserves this structure
    - every subset F of V_1 U E_1 U V_2 U E_2 is described via a dictionary that acts as an "indicator function"
"""

import networkx as nx
import numpy as np

from itertools import product
from collections import defaultdict

def f(F_description, wF = lambda x: 1):
    """
    Computes the function f(F) = sum_{v in F}w(v) + sum_{e in F}w(e).
    
    Arguments:
        - F_description: the description of F
        - wF: callable, a function that assigns a weight to every vertes/edge in F. Default: w(x) = 1 foreach x
    
    Returns:
        - the value of f(F)
    """
    return sum([wF(x) for x in F_description if F_description[x]])
    
def f_delta(g1, g2, intersect_m, m, u, v, wF = lambda x: 1):
    """
    Compute the increment in the value of f(g1 intersect(m') g2) where m' = m U {(u, v)}.
    See the description of f for further details.
                                                         
    It exploits the fact that the set g1 intersect(m') g2 can be computed by extending
    g1 intersect(m) g2 by adding:
        - u and v, if u and v have the same label
        - edges of the form u -> x from g1, if it is possible to find an edge v -> y with the same label in
            g2, with y in m(x)
            (NOTE: if the mapping m' is adding u -> x for the first time,
            it must be that the new edge has v as source, otherwise it would have been added by m)
        - edges of the form x -> u from g1, if it is possible to find an edge y -> v with the same label in
            g2, with y in m(x)
        - edges of the form v -> x from g2, if it is possible to find an edge u -> y with the same label in
            g1, with y in m(x)
        - edges of the form x -> v from g2, if it is possible to find an edge y -> u with the same label in
            g1, with y in m(x)
            
    At the same time, this function computes the extension of g1 intersect(m) g2 to g1 intersect(m') g2.
                                                                                                 
    Arguments:
        - g1, g2: the graphs to compare
        - intersect_m: the current description of g1 intersect(m) g2
        - m: the current relation m
        - u, v: the nodes to add to m to form m' = m U {(u, v)}
        - w: callable, a function that assigns a weight to every vertes/edge. Default: w(x) = 1 foreach x
        
    Returns:
        - delta: the increment f(g1 intersect(m') g2) - f(g1 intersect(m) g2)
        - intersection_delta: the description of (g1 intersect(m') g2) - (g1 intersect(m) g2)
    """
    delta = 0
    intersection_delta = {}
    
    # Add nodes
    if u[2] == v[2]: # u and v have the same label, they should be added to the intersection
                               # but first check whether they are not already there 
        delta += wF(u)*int(not intersect_m[u]) + wF(v)*int(not intersect_m[v])
        intersection_delta[u] = True
        intersection_delta[v] = True
        
    # Add edges
    
    for g_from, g_to, node_from, node_to in [(g1, g2, u, v), (g2, g1, v, u)]:
        
        V_from = g_from.nodes(data = "label")
        
        # u/v is the source
        for e in g_from.out_edges(node_from[1], data = "label"): # for each edge that has u/v as source (note: e = (source_id, target_id, label)) 
            
            if intersect_m[g_from, *e]: # if the edge is already in the intersecton, skip this edge, nothing to do
                continue
            
            for w in m[g_from, e[1], V_from[e[1]]]: # otherwise, for each mapping of the target of the edge in the other graph...
                if g_to.has_edge(node_to[1], w[1]): # if there is a corresponding edge that has v/u as source...
                    if g_to.get_edge_data(node_to[1], w[1])["label"] == e[2]: # and that edge has the same label as e
                        delta += wF(e) # then add its contribution...
                        
                        intersection_delta[g_from, *e] = True # and add it to the intersection
                        break
                    
        # u/v is the target (same as above with source/target roles reversed)
        for e in g_from.in_edges(node_from[1], data = "label"):
            
            if intersect_m[g_from, *e]:
                continue
            
            for w in m[g_from, e[0], V_from[e[0]]]:
                if g_to.has_edge(w[1], node_to[1]):
                    if g_to.get_edge_data(w[1], node_to[1])["label"] == e[2]:
                        delta += wF(e)
                        
                        intersection_delta[g_from, *e] = True # and add it to the intersection
                        break
    
    # print(intersection_delta)
    
    return delta, intersection_delta

def g(m, uv = None, wG = len):
    """
    Computes the function g(splits(m)) = sum_{v: |m(v)| > 1} w(v, m(v)).
    This function allows an optional parameter that allows to compute g(splits(m')), with
    m' = m U {(u,v)} without the need to update m.
    
    Arguments:
        - m: the relation m
        - uv: an optional tuple of nodes (u, v)
        - wG: callable, a function that assigns a weight to every split. Default: w(x) = len(x)
        
    Returns:
        - the value of g(splits(m))
    """
    
    if uv is not None:
        temp = {uv[0]: [uv[1]], uv[1]: [uv[0]]}
        uv = defaultdict(list, temp)
    else:
        uv = defaultdict(list)
        
    return sum([wG(m[x] + uv[x]) for x in m if len(m[x] + uv[x]) > 1]) # it is enough to iterate over the keys of m:
                                                                    #  -if neither u nor v are in m: these nodes are entirely new and it is the first
                                                                    #   time they are ever associated to anything, so they can't be split
                                                                    #  -if only u (equivalently, v) is in m, it means that v (u) was never associated
                                                                    #   before, hence only u can potentially be a split
                                                                    #
                                                                    # We also don't have to worry about duplicates in m[x] + uv[x]:
                                                                    #  suppose x = u, hence uv[u] = v. If m[u] already contains v, it
                                                                    #  means that (u,v) \in m, but then (u,v) could not have been considered
                                                                    #  as a candidate for extension.

def look_ahead(g1, g2, intersect_m_prime, u, v):
    """
    Computes the look_ahead(u, v) set.
    
    Let F = g1 intersect(m') g2, with m' = m U {(u,v)}
    This set is composed of:
        - all edges in g1 that are not in F, that have u as source/target, that have the same label as an edge in g2
            that has v as source/target
        - all edges in g2 that are not in F, that have v as source/target, that have the same label as an edge in g1
            that has u as source/target
            
    Arguments:
        - g1, g2: the graphs to compare
        - intersect_m_prime: the description of g1 intersect(m') g2
        - u, v: the nodes to add to m to form m' = m U {(u, v)}
    """
    out = defaultdict(bool)
    
    for g_from, g_to, node_from, node_to in [(g1, g2, u, v), (g2, g1, v, u)]:
        # in edges
        for e in g_from.in_edges(node_from, data = "label"): # for each edge that has u/v as source (Note: e = (source_id, target_id, label))
            
            if intersect_m_prime[g_from, *e]: # if the edge is already in the intersecton, skip this edge, nothing to do
                continue
            
            for e_cand in g_to.in_edges(node_to, data = "label"): # otherwise, for each edge that has v/u as source...
                if e_cand[2] == e[2]: # if the labels of the two edges coincide...
                    out[g_from, *e] = True # add the edge to the look_ahead set
                    break
                
        # out edges
        for e in g_from.out_edges(node_from, data = "label"): # for each edge that has u/v as source (Note: e = (source_id, target_id, label))
            
            if intersect_m_prime[g_from, *e]: # if the edge is already in the intersecton, skip this edge, nothing to do
                continue
            
            for e_cand in g_to.out_edges(node_to, data = "label"): # otherwise, for each edge that has v/u as source...
                if e_cand[2] == e[2]: # if the labels of the two edges coincide...
                    out[g_from, *e] = True # add the edge to the look_ahead set
                    break
                    
    return out

def similarity(g1, g2, wF = lambda x: 1, wG = len, verbose = False, random_state = 204475348541921399076734861309367813083):
    """
    Computes the similarity between g1 and g2.
    
    Arguments:
        - g1, g2: the two igraph.Graph objects to compare
        - wF: callable, weight function for f. Default: w(x) = 1 foreach x
        - wG: callable, weight function for g. Default: w(x) = len(x)
        - verbose: boolean, whether detailed output about the current state should be printed. Default: False
        - random_state: integer, a seed for the random number generator
        
    Returns:
        - the similarity score
        - the associated mapping
    """
    
    if(id(g1) == id(g2)):
        g2 = g2.copy()
    
    V1 = [(g1, *v) for v in g1.nodes(data = "label")]
    V2 = [(g2, *v) for v in g2.nodes(data = "label")]
    
    E1 = [(g1, *e) for e in g1.edges(data = "label")]
    E2 = [(g2, *e) for e in g2.edges(data = "label")]
    
    rng = np.random.default_rng(random_state)
    
    m = defaultdict(list)
    m_best = None
    
    intersect_m = defaultdict(bool) # dictionary with all the nodes and edges of g1 and g2:
                                    #  used to represent g1 intersect(m) g2
                                    #  and to update the score function with the minimal amount of computations
    
    curr_f = 0 # current value of f(intersect_m)
    best_score = 0
    
    candidates = list(product(V1, V2)) # set of nodes not yet included in m
    
    stop = False
    while not stop:
        
        # Compute candidates set
        candidates_loc = []
        max_increase = -float("inf")
        
        for u, v in candidates:
            delta_f, intersection_delta = f_delta(g1, g2, intersect_m, m, u, v, wF)
            loc_g = g(m, (u, v), wG)
            
            increase = delta_f - loc_g
            
            if increase == max_increase:
                candidates_loc += [((u,v), delta_f, loc_g, increase, intersection_delta)]
                
            if increase > max_increase:
                max_increase = increase
                candidates_loc = [((u,v), delta_f, loc_g, increase, intersection_delta)]
        
        candidates_loc_1 = []
        max_lookahead = -float("inf")
        stop = True
        for uv, delta_f, loc_g, increase, intersection_delta in candidates_loc:
            la = look_ahead(g1, g2, intersect_m | intersection_delta, uv[0], uv[1])
            fla = f(la, wF)
                
            if fla == max_lookahead:
                candidates_loc_1 += [(uv, delta_f, loc_g, increase, intersection_delta)]
                
                stop = stop and (fla <= 0) and (increase <= 0) # we can continue if lookahead is not empy or increase is positive 
                                                               # (since here we are only updating the candidates set,
                                                               # we also have to check that no other candidates allows for continuation)
                
            if fla > max_lookahead:
                max_lookahead = fla
                candidates_loc_1 = [(uv, delta_f, loc_g, increase, intersection_delta)]
                
                stop = (fla <= 0) and (increase <= 0) # we can continue if lookahead is not empy or increase is positive
                
        # Pick candidate
        chosen_candidate = rng.integers(0, len(candidates_loc_1))
        chosen_candidate = candidates_loc_1[chosen_candidate]
        
        uv, delta_f, loc_g, increase, intersection_delta = chosen_candidate
        u = uv[0]
        v = uv[1]
        
        # u = chosen_candidate[0][0]
        # v = chosen_candidate[0][1]
        # increase = chosen_candidate[1]
        # intersection_delta = chosen_candidate[2]
        
        # Update m/score(m)/g1 intersection(m) g2
        m[u] += [v]
        m[v] += [u]
        
        curr_f += delta_f
        curr_score = curr_f - loc_g
        
        if curr_score > best_score:
            best_score = curr_score
            m_best = { # deep copy of m
                    x: [y for y in mx] for x, mx in m.items()
                }
        
        intersect_m = intersect_m | intersection_delta
        # print(intersect_m)
        
        # Remove chosen candidates form pool of available candidates
        candidates.remove((u, v))
        
        if len(candidates) == 0:
            stop = True
        
        if verbose:
            print("\n------")
            print("Intersection:")
            print_intersection(intersect_m)
            print("..")
            print("m:")
            print_relation(m)
            print("Score:", curr_score, "f:", curr_f, "g:", loc_g, "increment:", increase, "stop:", stop, "id:", id(m), "id_best:", id(m_best))
    
    desc_g1_g2 = {
            **{v: True for v in list(V1) + list(V2)},
            **{e: True for e in list(E1) + list(E2)}
        }
    sim = best_score / f(desc_g1_g2, wF)
    
    return sim, m_best

def print_relation(m):
    graphs = set(id(v[0]) for v in m)
    
    for i, g in enumerate(graphs):
        print(f"From G_{i + 1} to G_{2 - i}:")
        
        for v, mv in m.items():
            if id(v[0]) == g:
                line = ", ".join([
                        f"({w[1]}, {w[2]})" for w in mv
                    ])
                
                print(f"\t({v[1]}, {v[2]}) -> {line}")

def print_intersection(intersect_m):
    graphs = set(id(x[0]) for x in intersect_m)

    for i, g in enumerate(graphs):
        print(f"G_{i + 1}:")
        
        line = []
        for x, included in intersect_m.items():
            if not included or id(x[0]) != g:
                continue
            
            if len(x) == 3: # node
                line += [f"({x[1]}, {x[2]})"]
            elif len(x) == 4:
                line += [f"({x[1]}, {x[2]}, {x[3]})"]
            
        print(", ".join(line))

if __name__ == "__main__":
    g1 = nx.DiGraph()
    g1.add_nodes_from([(i, {"label": str(i)}) for i in range(4)])
    g1.add_edges_from([(0, 1, {"label": "a"}), (1, 2, {"label": "b"}), (2, 3, {"label": "c"}), (1, 3, {"label": "d"})])
    
    g2 = nx.DiGraph()
    g2.add_nodes_from([(i, {"label": str(i)}) for i in range(5)])
    g2.add_edges_from([(0, 1, {"label": "a"}), (1, 2, {"label": "d"}), (2, 3, {"label": "c"}), (1, 3, {"label": "b"}), (4, 3, {"label": "e"})])
    
    sim, m = similarity(g1, g2, verbose = True)
    print("\n")
    print("Similarity:", sim)
    print("\n")
    print_relation(m)
