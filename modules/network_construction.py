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


# ===== Prepare =====
# taxa_koDict from make_table.py
taxa_koDict = None
with open(os.path.join(work_dir, "taxa_koDict.json"), 'r') as fp: taxa_koDict = json.load(fp)

# raw PLlist
PL_out = os.path.join(work_dir,"prelect_dir")
PLres = pd.read_csv(os.path.join(PL_out, "PLres.csv"), index_col=0)


# taxa_table
# taxa_table = pd.read_csv(os.path.join(work_dir, "taxa_table.csv"), index_col=0)
tmp_dir = os.path.join(work_dir, "tmp_files")


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

# ===== Main =====
edge = taxa2taxa(num_processes=threads)
edge.to_csv(os.path.join(tmp_dir, f"{target_label}_edge.csv"), index=False)
print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - network construction done!")
sys.exit(0)
