import pandas as pd
import numpy as np
import scipy
import warnings
from scipy.stats import fisher_exact
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=(SettingWithCopyWarning, FutureWarning))
import concurrent.futures
import argparse
import time
import json
import copy
import sys
import os

# input : work_dir, target_label(ex. 1_CRC), size_option(ex. degree, FC, KO)
# output : In Shiny, display graph and table (store in ./tmp_files first)

# parameters
work_dir = sys.argv[1] # taxa_table directory
target_label = sys.argv[2] # target label
p_val = float(sys.argv[3]) # p_val options
threads = int(sys.argv[4]) # threads
select_task = sys.argv[5] # selected classification task

# visualization refined options
size_opt = sys.argv[6] # size options
select_rank = sys.argv[7] # select rank
size_range_min = float(sys.argv[8]) # size range min
size_range_max = float(sys.argv[9]) # size range max
width_range_min = float(sys.argv[10]) # width range min
width_range_max = float(sys.argv[11]) # width range max




# ===== Prepare =====
script_dir = os.path.dirname(os.path.abspath(__file__))
# rank = ["phylum", "class", "order", "family", "genus", "none"]
sp_map = pd.read_csv(os.path.join(script_dir, "taxa_lineage", "species_map.csv")) # group info map
# taxa_koDict from make_table.py
taxa_koDict = None
with open(os.path.join(work_dir, "taxa_koDict.json"), 'r') as fp: taxa_koDict = json.load(fp)

# raw PLres
PL_out = os.path.join(work_dir,"prelect_dir")
PLres = pd.read_csv(os.path.join(PL_out, "PLres.csv"), index_col=0)

# taxa_table
tmp_dir = os.path.join(work_dir, "tmp_files")

# MC FC nessesary table
taxa_table = pd.read_csv(os.path.join(work_dir, "taxa_table.csv"), index_col=0)
meta = pd.read_csv(os.path.join(work_dir, "tmp_files", "meta.csv"), index_col=0)



# extract tendency taxa (tenedency PLres)
PLres_f = None
if select_task in ["BC", "Rg", "Cox"]: # binary classification, regression
    if target_label != "All":
        PLres_f = PLres[PLres["tendency"] == target_label].iloc[:, 0:3]
    else:
        PLres_f = PLres[PLres["selected"] == "Selected"].iloc[:, 0:3]

elif select_task == "MC": # multi classification
    tend_col = f"tendency_{target_label}"
    coef_col = f"coef_{target_label}"
    PLres_f = PLres[PLres[tend_col] == target_label].loc[:, ["FeatName", coef_col, tend_col]]


# summary
total_ko = []
for k in taxa_koDict:
    total_ko += taxa_koDict[k]
total_taxa = list(taxa_koDict.keys())
total_ko = list(set(total_ko))
PLres_f = PLres_f[PLres_f["FeatName"].isin(total_taxa)] # filter taxa without any ko

# collect tendency taxa total_ko background
feature_l = PLres_f["FeatName"].to_list()
ko_background = set()
for t in feature_l:
    ko_l = taxa_koDict[t]  # get KOs of taxa i
    for ko in ko_l:
        ko_background.add(ko)
ko_background = sorted(list(ko_background))


# ===== Main functions =====
def OR_module(matrix):
    """
    [[l1, r1], [l2, r2]]

    """
    # matrix
    L1, R1 = matrix[0][0], matrix[0][1]
    L2, R2 = matrix[1][0], matrix[1][1]

    if L1 == 0:
        dir_OR = 0
    elif (R1 > 0) and (L2 > 0):
        dir_OR = np.log10(L1/R1) - np.log10(L2/R2)
    else:
        # OR = np.inf if R1 == 0 else -(np.inf)
        # OR = np.inf
        dir_OR = 10 # if inf

    return dir_OR


# edge table
def calculate_edge(taxa_pair):
    global taxa_koDict, ko_background
    taxa1, taxa2 = taxa_pair
    mtable = pd.DataFrame({'KO': ko_background})  # make a dataframe, taxa1 taxa2 ko match
    mtable['taxa1'] = 0
    mtable['taxa1'][mtable['KO'].isin(taxa_koDict[taxa1])] = 1
    mtable['taxa2'] = 0
    mtable['taxa2'][mtable['KO'].isin(taxa_koDict[taxa2])] = 1

    # 2 by 2 table
    l1 = mtable[(mtable['taxa1'] == 1) & (mtable['taxa2'] == 1)].shape[0]
    l2 = mtable[(mtable['taxa1'] == 0) & (mtable['taxa2'] == 1)].shape[0]
    r1 = mtable[(mtable['taxa1'] == 1) & (mtable['taxa2'] == 0)].shape[0]
    r2 = mtable[(mtable['taxa1'] == 0) & (mtable['taxa2'] == 0)].shape[0]

    # fisher table
    m = [[l1, r1], [l2, r2]]
    res = fisher_exact(m, alternative='greater')
    weight = OR_module(m)
    return pd.DataFrame({'from': [taxa1], 'to': [taxa2], 'weight': [weight], 'p': [res.pvalue]})


# edge table
def taxa2taxa(num_processes):
    edgelist = []
    taxa_pairs = [(feature_l[i], feature_l[j]) for i in range(len(feature_l)) for j in range(i+1, len(feature_l))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        results = executor.map(calculate_edge, taxa_pairs)
    edgelist = list(results)

    edge = pd.concat(edgelist, axis=0, ignore_index=True)
    edge = edge[edge['weight'] != 0]
    # p.adjusted and filter
    edge["p_adj"] = scipy.stats.false_discovery_control(edge.p)
    edge = edge[edge["p_adj"] < p_val]

    # undirected attribute
    edge['Type'] = 'undirected'
    return edge



# Node table
unclassfied_species = []
def node_table(edge, opt):

    def ko_counts(taxa):
        ko_l = taxa_koDict[taxa]
        return len(ko_l)


    def fold_change(taxa, group):
        abund_list = taxa_table.loc[taxa,:]
        target_abund = abund_list[abund_list.index.isin(meta["Accession_ID"][meta["Labels"] == group])]
        other_abund = abund_list[abund_list.index.isin(meta["Accession_ID"][meta["Labels"] != group])]
        # FC
        fc = np.log2( (np.mean(target_abund)+1) / (np.mean(other_abund)+1) )
        return fc


    nodelist = []
    for n in np.unique(edge['from'].to_list() + edge['to'].to_list()):
        species_name = n.split('s__')[1]
        genus_name = n.split('.s__')[0].split('g__')[1]

        # Compatibility with selected task
        tendency = None
        if select_task in ["BC", "Rg", "Cox"]:
            tendency = PLres.loc[n, "tendency"]
        elif select_task == "MC":
            tendency = target_label

        # optional size
        if opt == "degree":
            size = edge[edge["from"] == n].shape[0] + edge[edge["to"] == n].shape[0] # degree of species
        elif opt == "FC":
            if select_task == "MC":
                fc = fold_change(n, target_label)
                abs_fc = abs(fc)
                size = abs_fc if abs_fc != np.inf else 10
            else:
                abs_fc = abs(PLres.loc[n,"logFC"])
                size = abs_fc if abs_fc != np.inf else 10
        elif opt == "KO":
            size = ko_counts(n)

        if select_rank in ["phylum", "class", "order", "family", "genus"]:
            subset = sp_map[select_rank][sp_map["genus"] == genus_name]
            if len(subset) == 0:
                global unclassfied_species
                unclassfied_species.append(n)
                edge = edge[~((edge['from'] == n) | (edge['to'] == n))] # remove node contains n
                continue
            else:
                group = group = subset.values[0]
            
            nodelist.append(pd.DataFrame({'id':[n],'label':[species_name],'sublabel':[genus_name], 'tendency':[tendency], 'raw_value':[size],'group':[group]}))
        else:
            nodelist.append(pd.DataFrame({'id':[n],'label':[species_name],'sublabel':[genus_name],'tendency':[tendency], 'raw_value':[size]}))

    node = pd.concat(nodelist, ignore_index=True)
    node = size_scaler(node)
    return node


def size_scaler(node_table):
    global size_range_min, size_range_max
    scale_node = copy.deepcopy(node_table)
    # Min-Max Scaler to scale size to 0-1
    min_value = scale_node['raw_value'].min()
    max_value = scale_node['raw_value'].max()
    scale_node['size'] = (scale_node['raw_value'] - min_value) / (max_value - min_value)

    # Scale to desired range (10-40)
    min_size = size_range_min
    max_size = size_range_max
    scale_node['size'] = min_size + scale_node['size'] * (max_size - min_size)
    return scale_node


def width_scaler(edge_table):
    global width_range_min, width_range_max
    scale_edge = copy.deepcopy(edge_table)
    # Min-Max Scaler to scale size to 0-1
    min_value = scale_edge['weight'].min()
    max_value = scale_edge['weight'].max()
    scale_edge['width'] = (scale_edge['weight'] - min_value) / (max_value - min_value)

    # Scale to desired range
    min_width = width_range_min
    max_width = width_range_max
    scale_edge['width'] = min_width + scale_edge['width'] * (max_width - min_width)
    return scale_edge


# ===== Main =====
edge = taxa2taxa(num_processes=threads)
# node table of edge
node = node_table(edge, size_opt)
# scale width of edge
scaling_edge = width_scaler(edge)

scaling_edge.to_csv(os.path.join(tmp_dir, f"{target_label}_edge.csv"), index=False)
node.to_csv(os.path.join(tmp_dir, f"{target_label}_node.csv"), index=False)
if len(unclassfied_species) != 0:
    with open(os.path.join(tmp_dir, f"{target_label}_unclassified_sp.txt"), "w") as fuc:
        for sp in unclassfied_species:
            fuc.write(f"{sp}\n")

print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - network construction done!")
sys.exit(0)
