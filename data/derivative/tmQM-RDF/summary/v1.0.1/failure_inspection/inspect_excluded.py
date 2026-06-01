# %% Locate root dir
from pydantic import constr
from sympy.printing.pretty.pretty_symbology import line_width
import os

ROOT_DIR = os.path.abspath(".")
while not ".prj_root" in os.listdir(ROOT_DIR):
    ROOT_DIR = os.path.abspath(os.path.join(ROOT_DIR, ".."))

# %% Header    
from matplotlib.ticker import MultipleLocator
from collections import Counter
from scipy import stats
from tqdm import tqdm

import matplotlib.transforms as pltran
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib
import time
import re

INPUT_FILES = {
		"viable": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v1.1", "viable_tmcs.txt"),
		"failures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v1.1", "tmQMg-L", "no_coverage.txt"),
		"ligands_atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v1.1", "tmQMg-L", "ligands_atoms_idx.csv"),
		"ligands_misc": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v74k", "ligands_misc_info.csv"),
		"ligands_fingerprints": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v74k", "ligands_fingerprints.csv"),
		"ligands_descriptors": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg-L", "v74k", "ligands_descriptors.csv"),
		"periodic_table": os.path.join(ROOT_DIR, "data", "raw", "pubChem", "data", "PubChemElements_all.csv"),
		"tmQM_series_properties": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "misc", "tmQM_series_property_table.csv"),
        "tmQM_series_property_descriptions": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "misc", "tmQM_series_property_descriptions.csv"), # TODO
		"tmQM_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v1.1", "tmQM", "uNatQ_graphs"),
		"tmQMg_graphs": os.path.join(ROOT_DIR, "data", "raw", "tmQMseries", "data", "tmQMg", "v74.548k", "u-NatQ_graphs"),
		"tmQMg_L_graphs": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "intermediate", "v1.1", "tmQMg-L", "uNatQ_graphs"),
        "details": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF-1Ksel", "intermediate", "v1.1", "details.csv")
	}

OUTPUT_FILES = {
		"centres": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "tmQMg_centres.txt"),
		"coverage_a": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "coverage_a.txt"),
		"coverage_l": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "coverage_l.txt"),
        "atoms": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "atoms.txt"),
        "latmods": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "latmods.txt"),
        "electronic": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "electronic.txt"),
		"figures": os.path.join(ROOT_DIR, "data", "derivative", "tmQM-RDF", "summary", "v1.1", "failure_inspection", "figures")
	}

INSPECT_CENTRES = False
INSPECT_COVERAGE = False
INSPECT_ATOMS = True
INSPECT_LIGAND_ATTACHMENT = False
INSPECT_TMC_PROPERTIES = False

def checkpoint(filename):
    """
    Decorator factory that produces a decorator with the following effect:
        - if the file 'filename' exists, then the decorated function is NOT invoked and the content of 'filename',
            evalueted via eval(), is returned instead
        - otherwise, the decorated function is invoked and its output is stored in 'filename' before being returned
    """

    def checkpoint(func):

        def decorated_function():
            if os.path.exists(filename):
                with open(filename, "r") as f:
                    return eval(f.read())
                
            out = func()

            with open(filename, "w") as f:
                f.write(str(out))

            return out
        
        return decorated_function
    
    return checkpoint

# %% Centre distribution ----

@checkpoint(OUTPUT_FILES["centres"])
def retrieve_centres_from_failures():
    tmqmg = [x for x in os.listdir(INPUT_FILES["tmQMg_graphs"]) if x.endswith(".gml")]

    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")

    centres = []
    centres_failures = []

    for tmc in tqdm(tmqmg):
        with open(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc), "r") as g:
            centre = g.read().split("\n")[6].strip().split(" ")[1][1:-1]

        centres.append(centre)

        if tmc.split(".")[0] in failures:
            centres_failures.append(centre)

    return {"centres": centres, "centres_failures": centres_failures}

def inspect_metal_centres(centres_dict):
    centres = centres_dict["centres"]
    centres_failures = centres_dict["centres_failures"]

    pubchem = pd.read_csv(INPUT_FILES["periodic_table"], sep = ",")

    # Centres distribution (normalised) in tmQMg
    count = Counter(centres)
    count = {
        el: c 
        for el, c in sorted(count.items(), key = lambda ec: pubchem.loc[pubchem["Symbol"] == ec[0], "AtomicNumber"].iloc[0])
    }
    Z = sum(count.values())
    count = {el: c/Z for el, c in count.items()}

    # Centres distribution (normalised) in excluded TMCs
    count_fail = Counter(centres_failures)
    count_fail = {
        el: c 
        for el, c in sorted(count_fail.items(), key = lambda ec: pubchem.loc[pubchem["Symbol"] == ec[0], "AtomicNumber"].iloc[0])
    }
    Z = sum(count_fail.values())
    count_fail = {el: c/Z for el, c in count_fail.items()}
    
    # Centres distribution (normalised) in tmQM-RDF
    sel_info = pd.read_csv(INPUT_FILES["details"], index_col = 0)
    
    
    count_ds = Counter(sel_info.loc[:, "core"])
    count_ds = {
        el: c 
        for el, c in sorted(count_ds.items(), key = lambda ec: pubchem.loc[pubchem["Symbol"] == ec[0], "AtomicNumber"].iloc[0])
    }
    Z = sum(count_ds.values())
    count_ds = {el: c/Z for el, c in count_ds.items()}

    # Ensure both count dicts cover the same elements
    for el in count_fail:
        if el not in count_ds:
            count_ds[el] = 0


    for el in count_ds:
        if el not in count_fail:
            count_fail[el] = 0

    # Ratios
    ratio = {
        el: count_fail[el]/c for el, c in count_ds.items()
    }

    # Plot normalised counts
    fig, axs = plt.subplots(3, 1, figsize = (6.2, 4.8))
    x = np.arange(len(count_fail))
    
    width = 0.7

    for ax, cnt, col, nm in zip(axs, [count, count_fail, count_ds], ["green", "red", "blue"], ["tmQMg", "failures", "tmQM-RDF"]):
        ## Plot tmQMg
        ax.bar(x, cnt.values(), width = width, color = col, alpha = 0.3)

        ## Plot failures
        #ax.bar(x, count_fail.values(), width = width, color = "red", alpha = 0.3)

        ## Plot tmQM-RDF
        #ax.bar(x, count_ds.values(), width = width, color = "blue", alpha = 0.3)

        ## Setup xticks
        ax.set_xticks(x, count_fail.keys(), fontsize = "x-small")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
        
        dx = -2/72.
        dy = 3/72. # shifts (in inches)
        offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        
        for label in ax.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)
        
        ## Axis labels
        ax.set_xlabel("Metal centres (ordered by atomic number)")
        ax.set_ylabel(f"Normalised counts ({nm})")

    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "centres.pdf"), format = "pdf", bbox_inches = "tight")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "centres.png"), format = "png", bbox_inches = "tight", dpi = 200)

    # Plot normalised counts as color bands
    fig, axs = plt.subplots(3, 1, figsize = (6.2, 4.8))
    x = np.arange(len(count_fail))
    
    width = 0.7

    for ax, cnt, nm in zip(axs, [count, count_fail, count_ds], ["tmQMg", "failures", "tmQM-RDF"]):
        ## Plot tmQMg
        #ax.bar(x, cnt.values(), width = width, color = col, alpha = 0.3)
        img = ax.imshow(np.array(list(cnt.values())).reshape(1, -1), cmap = "magma", aspect = "auto")

        ## Plot failures
        #ax.bar(x, count_fail.values(), width = width, color = "red", alpha = 0.3)

        ## Plot tmQM-RDF
        #ax.bar(x, count_ds.values(), width = width, color = "blue", alpha = 0.3)

        ## Setup xticks
        ax.set_xticks(x, count_fail.keys(), fontsize = "x-small")
        plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
        
        dx = -2/72.
        dy = 3/72. # shifts (in inches)
        offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        
        for label in ax.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)
        
        ## Axis labels
        ax.set_xlabel("Metal centres (ordered by atomic number)")
        ax.set_ylabel(nm)

        ax.set_yticks([])

        cbar = plt.colorbar(img, ax = ax, orientation = "vertical")

    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "dists_cbar.pdf"), format = "pdf", bbox_inches = "tight")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "dists_cbar.png"), format = "png", bbox_inches = "tight", dpi = 200)

    # Plot ratios
    fig, ax = plt.subplots(figsize = (6.2, 4.8))
    x = np.arange(len(count_fail))
    
    ax.bar(x, ratio.values(), color = "orange")

    ax.hlines(1, min(x), max(x), colors = "orange", linestyles = "dashed")
    #lb = min(rel_diff.values())

    #ax.vlines(x, lb*1.1, 0, colors = "orange", linestyles = "dashed")

    ## Setup xticks
    ax.set_xticks(x, count_fail.keys(), fontsize = "x-small")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
    
    dx = -2/72.
    dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    
    ## Axis labels
    ax.set_xlabel("Metal centres (ordered by atomic number)")
    ax.set_ylabel("Ratio of normalised counts (fails/tmQM-RDF)")

    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "ratios.pdf"), format = "pdf", bbox_inches = "tight")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "ratios.png"), format = "png", bbox_inches = "tight", dpi = 200)

# %% Ligand coverage inspection ----

#@checkpoint(OUTPUT_FILES["coverage_l"])
def retireve_ligand_coverage_ligands():
 
    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")
    
    subgraph_data = pd.read_csv(INPUT_FILES["ligands_atoms"], sep = ",", header = 0)
    c = Counter([x.split("-")[0] for x in subgraph_data.iloc[:,0]])

    tmcs = []

    for tmc in tqdm(failures):
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], f"{tmc}.gml"))
        g.remove_nodes_from([n for n, attr in g.nodes(data = "node_label") if attr == g.graph["meta_data"]["metal_center_element"]])

        tmQMg_nl = nx.number_connected_components(g)
        tmQMgL_nl = c[tmc]

        tmcs += [(tmQMg_nl, tmQMgL_nl)]

        if tmQMg_nl == tmQMgL_nl:
            print(tmc)

    return tmcs

@checkpoint(OUTPUT_FILES["coverage_a"])
def retireve_ligand_coverage_atoms():
    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")
    
    subgraph_data = pd.read_csv(INPUT_FILES["ligands_atoms"], sep = ",", header = 0)
    subgraph_data["subgraph"] = subgraph_data["subgraph"].apply(lambda x: x.split("-")[0])

    tmcs = []

    for tmc in tqdm(failures):
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], f"{tmc}.gml"))
        g.remove_nodes_from([n for n, attr in g.nodes(data = "node_label") if attr == g.graph["meta_data"]["metal_center_element"]])

        tmQMg_nat = len(g)
        tmQMgL_nat = 0

        tmc_data = subgraph_data.loc[subgraph_data["subgraph"] == tmc, :]
        for i in range(tmc_data.shape[0]):
            tmQMgL_nat += len(tmc_data.iloc[i, 1].strip().split(" "))

        tmcs += [(tmQMg_nat, tmQMgL_nat)]

    return tmcs


def plot_ligand_coverage(coverage_a, coverage_l):
    #fig, ax = plt.subplots(figsize = (4.8, 4.8))
    
    for coverage, labels, ml in zip([coverage_a, coverage_l], ["non-metal atoms", "ligands"], [20, 2]):
        fig, ax = plt.subplots()

        tmQMgL = [c[1] for c in coverage]
        tmQMg = [c[0] for c in coverage]

        ax.scatter(tmQMg, tmQMgL, color = "orange", alpha = 0.5, s = 3)
        ax.plot([0, max(tmQMg + tmQMgL)], [0, max(tmQMg + tmQMgL)], color = "grey", linewidth = 0.5, linestyle = "dashed")

        ax.set_xlabel(f"No. of {labels} in tmQMg")
        ax.set_ylabel(f"No. of {labels} in tmQMg-L")

        ax.xaxis.set_major_locator(MultipleLocator(ml))
        ax.yaxis.set_major_locator(MultipleLocator(ml))

        ax.set_aspect("equal")

        key = labels.split(" ")[-1]
        fig.savefig(os.path.join(OUTPUT_FILES["figures"], f"coverage_{key}.pdf"), format = "pdf", bbox_inches = "tight")
        fig.savefig(os.path.join(OUTPUT_FILES["figures"], f"coverage_{key}.png"), format = "png", bbox_inches = "tight", dpi = 200)

# %% Atoms and elements ----
@checkpoint(OUTPUT_FILES["atoms"])
def retrieve_atoms_and_sizes():
    tmqmg = [x for x in os.listdir(INPUT_FILES["tmQMg_graphs"]) if x.endswith(".gml")]

    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")

    out = {
        "atoms": {
            "tmQMg": dict(),
            "failures": dict(),
            "tmQM-RDF": dict()
        },
        "sizes": {
            "tmQMg": list(),
            "failures": list(),
            "tmQM-RDF": list()
        }
    }

    for tmc in tqdm(tmqmg):
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc))
        s = g.graph["meta_data"]["n_atoms"]
        c = g.graph["meta_data"]["element_counts"]

        key = "failures" if tmc.split(".")[0] in failures else "tmQM-RDF"

        out["sizes"]["tmQMg"].append(s)
        out["sizes"][key].append(s)

        for el, count in c.items():
            out["atoms"]["tmQMg"][el] = out["atoms"]["tmQMg"].get(el, 0) + int(count)
            out["atoms"][key][el] = out["atoms"][key].get(el, 0) + int(count)

    return out

def plot_atoms_and_sizes(atoms_and_sizes, metals = []):
    pubchem = pd.read_csv(INPUT_FILES["periodic_table"], sep = ",")

    atoms = atoms_and_sizes["atoms"]
    sizes = atoms_and_sizes["sizes"]

    nice_labels = {
        "tmQMg": "tmQMg",
        "failures": "Excluded TMCs",
        "tmQM-RDF": "tmQM-RDF"
    }

    keys = ["tmQMg", "failures", "tmQM-RDF"]
    for k in keys:
        for a in atoms[k]:
            for j in keys:
                if j == k:
                    continue
                if a not in atoms[j]:
                    atoms[j][a] = 0
        
    for a in atoms["tmQMg"]:
        if a not in atoms["failures"]:
            atoms["failures"][a] = 0

    # Plot element counts
    
    """Sort elements by atomic number (metal centres first, then all the rest)"""
    def sort_by_atnum(el_list):
        return sorted(el_list, key = lambda e: pubchem.loc[pubchem["Symbol"] == e, "AtomicNumber"].iloc[0])

    sorted_metals = sort_by_atnum(metals)
    sorted_nonmetals = sort_by_atnum([e for e in atoms["tmQMg"] if e not in metals])
    sorted_elements = sorted_metals + sorted_nonmetals

    fig, ax = plt.subplots(figsize = (12.4, 4.8), layout = "constrained")
    x = np.arange(len(atoms["tmQMg"]))*2
    x[len(sorted_metals):] += 1

    width = 0.5
    for mult, key, col, hatchpar, in zip([-1, 0, 1], keys, ["white", "red", "blue"], [{"edgecolor": "black", "hatch": "//"}, {}, {}]):
        # Centres distribution (normalised)
        count = {
            el: atoms[key][el] 
            for el in sorted_elements
        }
        Z = sum(count.values())
        
        # Plot bars
        ax.bar(x + width*mult, np.array(list(count.values()))/Z, width = width, color = col, linewidth = 0.5, **hatchpar, label = nice_labels[key])
    
    ax.set_yscale("log")
    ax.legend(loc = "upper left")

    # Add metal/non-metal separator
    boundary_idx = len(sorted_metals)
    boundary_x = 0.5*(x[boundary_idx] + x[boundary_idx - 1])
    ax.vlines(x = boundary_x, ymin = 0, ymax = 1, transform = ax.get_xaxis_transform(), linewidth = 1, color = "k")

    ax.text(boundary_x - 0.5, 0.96, "Metals", transform = ax.get_xaxis_transform(), ha = "right")
    ax.text(boundary_x + 0.5, 0.96, "Non metals", transform = ax.get_xaxis_transform(), ha = "left")

    ## Setup xticks
    ax.set_xticks(x, count.keys(), fontsize = "x-small")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
    
    dx = -2/72.
    dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    
    ## Axis labels
    ax.set_xlabel("Elements (ordered by atomic number)")
    ax.set_ylabel("Frequencies (log scale)")

    #plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "elements.pdf"), format = "pdf")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "elements.png"), format = "png", dpi = 200)

    # Plot sizes
    fig, axs = plt.subplots(1, 3, figsize = (12.4, 4.8), layout = "constrained", sharey = True)
    
    for ax, key, col, hathchpar in zip(axs, keys, ["white", "red", "blue"], [{"edgecolor": "black", "hatch": "//"}, {}, {}]):
        # Centres distribution (normalised) in tmQMg
        
        ## Plot tmQMg
        ax.hist(sizes[key], color = col, density = True, **hathchpar, label = nice_labels[key])
        ax.legend(loc = "upper left")

        ## Setup xticks
        #ax.set_xticks(x, count.keys(), fontsize = "x-small")
        #plt.setp(ax.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
        
        #dx = -2/72.
        #dy = 3/72. # shifts (in inches)
        #offset = pltran.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        
        #for label in ax.xaxis.get_majorticklabels():
        #    label.set_transform(label.get_transform() + offset)
        
        ## Axis labels
        ax.set_xlabel("TMC size (no. of atoms)")
    
    axs[0].set_ylabel("Density")

    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "sizes.pdf"), format = "pdf")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "sizes.png"), format = "png", dpi = 200)


# %% Ligand attachments ----
def compute_latmod(g):
    # Find metal centre and its closure
    mc = [n for n, attr in g.nodes(data = "node_label") if attr == g.graph["meta_data"]["metal_center_element"]][0]
    closure = set(list(g[mc].keys()) + [mc])

    # Identify ligands (as connected components of 'g - mc')
    g_ligs = g.copy()
    g_ligs.remove_node(mc)
    ligs = dict(enumerate(nx.connected_components(g_ligs)))
    
    # identify ligand attachment modes
    latmod = dict()
    for lig, nodes in ligs.items():
        # Find binding atoms (nodes in ligand that also belong to the closure of the mc)
        binding_nodes = g.subgraph(closure & nodes)

        # Find all connected components among the binding atoms inside ligand
        binding_sites = list(nx.connected_components(binding_nodes))

        # Relabel binding atoms from nx id to chemical label
        binding_sites = [[g.nodes[n]['node_label'] for n in site] for site in binding_sites]

        # store binding sites: this is a list of lists
        # where each inner list is a set of atoms that are adiacent and bind to the mc
        # (it can happen that a binding site list is empty due to problems in tmQMg encoding (isolated nodes in the graph))
        if len(binding_sites) > 0:
            latmod[lig] = binding_sites

    return latmod

@checkpoint(OUTPUT_FILES["latmods"])
def retrieve_ligand_attachments():
    """
    In this function, ligand attachments are retrieved as strings of the form
    x.y.z. ...
    where x, y, z, ... are integers. Dots separate positions. The integer in position n
    counts how many binding sites of size n are present.

    E.g.: if the latmod of a ligand looks is [[C], [C,C,C], [C, C, C]], it will
        be encoded as 1.0.2
    """
    def encode(latmod):
        counts = Counter(len(site) for site in latmod)

        max_size = max(counts.keys())
        padded_counts = [
            counts[i] if i in counts else 0 
            for i in range(1, max_size + 1)
        ]

        return ".".join([str(c) for c in padded_counts])

    tmqmg = [x for x in os.listdir(INPUT_FILES["tmQMg_graphs"]) if x.endswith(".gml")]

    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")

    latmods = {
        "tmQMg": [],
        "failures": [],
        "tmQM-RDF": []
    }

    for tmc in tqdm(tmqmg):
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc))
    
        gname = tmc.split(".")[0]
        key = "failures" if gname in failures else "tmQM-RDF"

        loc_latmods = compute_latmod(g)

        for lm in loc_latmods.values():
            elm = encode(lm)
    
            latmods["tmQMg"].append(elm)
            latmods[key].append(elm)

    return latmods

def plot_latmods(latmods):
    full_token_set = set(sum(list(latmods.values()), []))

    counters = {
        k: Counter(v) for k, v in latmods.items()
    }

    for k in ["tmQMg", "tmQM-RDF", "failures"]:
        for tok in full_token_set:
            if tok not in counters[k]:
                counters[k][tok] = 0

    sorted_tokens = sorted(list(full_token_set))

    # Plot
    fig, axs = plt.subplots(3, 1, figsize = (6.2, 4.8))
    x = np.arange(len(sorted_tokens))
    
    width = 0.7

    for ax, key, col in zip(axs, ["tmQMg", "failures", "tmQM-RDF"], ["green", "red", "blue"]):
        y = [counters[key][tok] for tok in sorted_tokens]
        ax.bar(x, y, width = width, color = col, alpha = 0.3)
        ax.set_yscale("log")

        ## Plot failures
        #ax.bar(x, count_fail.values(), width = width, color = "red", alpha = 0.3)

        ## Plot tmQM-RDF
        #ax.bar(x, count_ds.values(), width = width, color = "blue", alpha = 0.3)

        ## Setup xticks
        ax.set_xticks([])
        ## Axis labels
        ax.set_xlabel("Binding modalities (lexycographic order)")
        ax.set_ylabel(f"Counts ({key})")

    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "latmods.pdf"), format = "pdf", bbox_inches = "tight")
    fig.savefig(os.path.join(OUTPUT_FILES["figures"], "latmods.png"), format = "png", bbox_inches = "tight", dpi = 200)

# %% Electronic properties ----
@checkpoint(OUTPUT_FILES["electronic"])
def retrieve_electronic_properties():
    tmqmg = [x for x in os.listdir(INPUT_FILES["tmQMg_graphs"]) if x.endswith(".gml")]

    with open(INPUT_FILES["failures"], "r") as f:
        failures = f.read().split("\n")

    # Retrieve property names
    blacklist = ["meta_data", "feature_n_atoms"]

    temp = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmqmg[0]))
    pnames = [n for n in temp.graph.keys() if n not in blacklist]

    properties = {
        "tmQMg": {p: [] for p in pnames},
        "failures": {p: [] for p in pnames},
        "tmQM-RDF": {p: [] for p in pnames}
    }

    for tmc in tqdm(tmqmg):
        g = nx.read_gml(os.path.join(INPUT_FILES["tmQMg_graphs"], tmc))

        gname = tmc.split(".")[0]
        key = "failures" if gname in failures else "tmQM-RDF"

        for p in pnames:
            properties["tmQMg"][p].append(g.graph[p])
            properties[key][p].append(g.graph[p])

    return properties

def check_distributions(properties):
    discrete = ["feature_n_electrons", "feature_charge"]
    continuous = [p for p in properties["tmQMg"] if p not in discrete]

    print("KS test on continuous electronic properties:")
    width = 0.5
    x0 = np.arange(len(continuous))*2
    fig0, ax0 = plt.subplots(figsize = (12.4, 4.8), layout = "constrained")
    for i, p in enumerate(continuous):
        x = properties["tmQM-RDF"][p]
        y = properties["failures"][p]
        z = properties["tmQMg"][p]
        """

        test = stats.kstest(x, y)

        print(f"\tProperty: {p} || tmQM-RDF/failures: {np.mean(x)} / {np.mean(y)} || p-value = {test.pvalue:.6f}")
        fig, ax = plt.subplots()
        ax.hist(y, 100, color = "red", alpha = 0.6, label = "failures", density = True)
        ax.hist(x, 100, color = "blue", alpha = 0.6, label = "tmQM-RDF", density = True)
        ax.set_xlabel(p)
        ax.legend()
        
        fig.savefig(os.path.join(OUTPUT_FILES["figures"], "properties", f"{p}.pdf"), format = "pdf", bbox_inches = "tight")
        fig.savefig(os.path.join(OUTPUT_FILES["figures"], "properties",f"{p}.png"), format = "png", bbox_inches = "tight", dpi = 200)
        """

        s = np.std(z)
        ax0.bar(x0[i] - width, [np.mean(z/s)], color = "white", width = width, hatch = "//", edgecolor = "black", yerr = np.std(z/s), label = "tmQMg")
        ax0.bar(x0[i], [np.mean(y/s)], color = "red", width = width, yerr = np.std(y/s), label = "Excluded TMCs")
        ax0.bar(x0[i] + width, [np.mean(x/s)], color = "blue", width = width, yerr = np.std(x/s), label = "tmQM-RDF")

    ax0.hlines(0, 0, 1, transform = ax0.get_yaxis_transform(), color = "black")
    
    ## Setup xticks
    ax0.set_xticks(
        x0, 
        [l.replace("feature_", "").replace("target_", "").replace("tzvp_", "") for l in continuous], 
        fontsize = "x-small"
    )
    plt.setp(ax0.xaxis.get_majorticklabels(), rotation = 60, rotation_mode = "anchor", ha = "right")
    
    dx = -2/72.
    dy = 3/72. # shifts (in inches)
    offset = pltran.ScaledTranslation(dx, dy, fig0.dpi_scale_trans)
    
    for label in ax0.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    
    ## Axis labels
    ax0.set_xlabel("Electronic properties")
    ax0.set_ylabel("Values (scaled)")

    # Handle legend
    handles, labels = ax0.get_legend_handles_labels()
    ax0.legend(handles[:3], labels[:3], loc = "upper left")

    plt.tight_layout()
    fig0.savefig(os.path.join(OUTPUT_FILES["figures"], "property_summary.pdf"), format = "pdf")
    fig0.savefig(os.path.join(OUTPUT_FILES["figures"], "property_summary.png"), format = "png", dpi = 200)

    for p in discrete:
        # Plot
        fig, axs = plt.subplots(3, 1, figsize = (6.2, 4.8))
        
        width = 0.7

        for ax, key, col in zip(axs, ["tmQMg", "failures", "tmQM-RDF"], ["green", "red", "blue"]):
            ax.hist(properties[key][p], width = width, color = col, alpha = 0.3)
            # ax.set_yscale("log")

            ## Plot failures
            #ax.bar(x, count_fail.values(), width = width, color = "red", alpha = 0.3)

            ## Plot tmQM-RDF
            #ax.bar(x, count_ds.values(), width = width, color = "blue", alpha = 0.3)

            ## Axis labels
            ax.set_xlabel(p)
            ax.set_ylabel(f"({key})")

        fig.savefig(os.path.join(OUTPUT_FILES["figures"], "properties", f"{p}.pdf"), format = "pdf", bbox_inches = "tight")
        fig.savefig(os.path.join(OUTPUT_FILES["figures"], "properties",f"{p}.png"), format = "png", bbox_inches = "tight", dpi = 200)

# %% Main statement ----
if __name__ == "__main__":
    # Centre distribution
    centres_dict = None
    if INSPECT_CENTRES:
        centres_dict = retrieve_centres_from_failures()
        inspect_metal_centres(centres_dict)

    if INSPECT_COVERAGE:
        coverage_a = retireve_ligand_coverage_atoms()
        coverage_l = retireve_ligand_coverage_ligands()
        plot_ligand_coverage(coverage_a, coverage_l)

    if INSPECT_ATOMS:
        if centres_dict is None:
            centres_dict = retrieve_centres_from_failures()

        atoms_and_sizes = retrieve_atoms_and_sizes()
        plot_atoms_and_sizes(atoms_and_sizes, np.unique(centres_dict["centres"]))

    if INSPECT_LIGAND_ATTACHMENT:
        latmods = retrieve_ligand_attachments()
        plot_latmods(latmods)
    
    if INSPECT_TMC_PROPERTIES:
        properties = retrieve_electronic_properties()
        check_distributions(properties)