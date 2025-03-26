#!/usr/bin/env python

import pandas as pd
import numpy as np
import scipy
import warnings
from scipy.stats import fisher_exact
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=(SettingWithCopyWarning, FutureWarning))
import concurrent.futures
# import pyreadr
import argparse
import time
import json
import copy
import sys
import os

# input : work_dir, target_label(ex. 1_CRC), p_val
# ouput : In Shiny, display vital pathway table (store in ./tmp_files first)

# parameters
work_dir = sys.argv[1] # taxa_table directory
target_label = sys.argv[2] # target label
p_val = float(sys.argv[3]) # p_val options
threads = int(sys.argv[4]) # threads

# ===== Prepare =====
script_dir = os.path.dirname(os.path.abspath(__file__))
tmp_dir = os.path.join(work_dir, "tmp_files")
# map_id to Kos
# pathwayInfo from ko_search module
whole_komap_dir = os.path.join(script_dir, "whole_komap", "total_KO_pathway.txt")
pathwayInfo = pd.read_csv(whole_komap_dir, sep="\t")
# whole_komap_dir = os.path.join(script_dir, "modules", "whole_komap", "total_KO_pathway.rds")
# pathwayInfo = pyreadr.read_r(whole_komap_dir)
# pathwayInfo = pathwayInfo[None]


# process pathway dict
taxa_koDict = None
with open(os.path.join(work_dir, "taxa_koDict.json"), 'r') as fp: taxa_koDict = json.load(fp) # taxa_koDict from make_table.py

path_koDict = {} # { mapid : [ko1, ko2, ..]}
for i in range(1, pathwayInfo.shape[0]+1):
    if pathwayInfo['mapid'][i] in path_koDict:
        path_koDict[pathwayInfo['mapid'][i]].append(pathwayInfo['KO'][i])
    else:
        path_koDict[pathwayInfo['mapid'][i]] = [pathwayInfo['KO'][i]]
for i in path_koDict:
    path_koDict[i] = list(set(path_koDict[i]))

# clean KO background (taxa_ko, path_ko)
taxa_table = pd.read_csv(os.path.join(work_dir, "taxa_table.csv"), sep=",", index_col=0)
# collect tendency taxa total_ko background
feature_l = set(taxa_table.index.to_list()).intersection(set(taxa_koDict.keys()))
ko_background = set()
for t in feature_l:
    ko_l = taxa_koDict[t]  # get KOs of taxa i
    for ko in ko_l:
        ko_background.add(ko)
ko_background = sorted(list(ko_background))
# taxa_koDict
filtered_taxa_koDict = {}
for k1, v1 in taxa_koDict.items():
    filterd_val = [ko for ko in v1 if ko in ko_background]
    filtered_taxa_koDict[k1] = filterd_val
taxa_koDict = filtered_taxa_koDict # update
# path_koDict
filtered_path_koDict = {}
for k2, v2 in path_koDict.items():
    filterd_val = [ko for ko in v2 if ko in ko_background]
    filtered_path_koDict[k2] = filterd_val
path_koDict = filtered_path_koDict # update


# get pathway taxa list
# when taxa_ko and path_ko are intersection and ko is known operon...
op_ko = None
with open(os.path.join(script_dir, "whole_komap", "op_ko.txt"), "r") as op_f:
    op_ko = [line for line in op_f.read().splitlines()]

path_taxaDict = {} # restore
for pk, pv in path_koDict.items():
    path_taxaDict[pk] = []
    # browse taxa_koDict
    for tk, tv in taxa_koDict.items():
        if ( set(pv).intersection(set(tv)) ) and ( set(tv).intersection(set(op_ko)) ):
            path_taxaDict[pk].append(tk)

# summary
total_taxa = list(taxa_koDict.keys())
total_ko = []
for k in taxa_koDict:
    total_ko += taxa_koDict[k]
total_ko = list(set(total_ko))


# ===== Main functions =====
def calculate_edge(taxa_pair):
    taxa1, taxa2 = taxa_pair
    mtable = pd.DataFrame({'ID': total_ko})  # make a dataframe, taxa1 taxa2 pairwise ko match
    mtable['taxa1'] = 0
    mtable['taxa1'][mtable['ID'].isin(taxa_koDict[taxa1])] = 1
    mtable['taxa2'] = 0
    mtable['taxa2'][mtable['ID'].isin(taxa_koDict[taxa2])] = 1

    # fisher table
    l1 = mtable[(mtable['taxa1'] == 1) & (mtable['taxa2'] == 1)].shape[0]
    r1 = mtable[(mtable['taxa1'] == 0) & (mtable['taxa2'] == 1)].shape[0]
    l2 = mtable[(mtable['taxa1'] == 1) & (mtable['taxa2'] == 0)].shape[0]
    r2 = mtable[(mtable['taxa1'] == 0) & (mtable['taxa2'] == 0)].shape[0]
    m = [[l1, r1], [l2, r2]]
    res = fisher_exact(m, alternative='greater')
    return pd.DataFrame({'from': [taxa1], 'to': [taxa2], 'weight': [res.statistic], 'p': [res.pvalue]})


def make_whole_edge(num_processes):
    whole_edgelist = []
    taxa_pairs = [(total_taxa[i], total_taxa[j]) for i in range(len(total_taxa)) for j in range(i+1, len(total_taxa))]

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        results = executor.map(calculate_edge, taxa_pairs)
        for result in results:
            whole_edgelist.append(result)

    # clean table
    whole_edge = pd.concat(whole_edgelist, axis=0, ignore_index=True)
    whole_edge = whole_edge[np.isfinite(whole_edge['weight'])]
    whole_edge = whole_edge[whole_edge['weight'] > 0]
    whole_edge['Type'] = 'undirected'

    # p.adjusted and filter
    whole_edge['p_adj'] = scipy.stats.false_discovery_control(whole_edge.p)
    whole_edge = whole_edge[whole_edge['p_adj'] < p_val]
    whole_edge['ID'] = list(range(whole_edge.shape[0]))  # assign ID
    whole_edge.to_csv(os.path.join(NC_dir, "whole_taxa2taxa_edge.tsv"), sep="\t", index=False)
    return


def vital_pathway(filter_edge, filter_node, backgrond_edge):
    # ===== Step 1 =====
    network_taxa = filter_node["id"].to_list()
    not_network_taxa = [t for t in total_taxa if t not in network_taxa]

    # browse each pathway & calculate
    calres_list = []
    for path, path_taxa_list in path_taxaDict.items():
        # 非 pathway i 的 taxa
        taxa_not_in_list = [t for t in total_taxa if t not in path_taxa_list]
        # value
        l1 = len(set(network_taxa).intersection(path_taxa_list))
        l2 = len(set(network_taxa).intersection(taxa_not_in_list))
        r1 = len(set(not_network_taxa).intersection(path_taxa_list))
        r2 = len(set(not_network_taxa).intersection(taxa_not_in_list))
        m = [[l1, r1],[l2, r2]]
        res = fisher_exact(m, alternative='greater')
        calres_list.append(pd.DataFrame({'pathway':[path], 'step1_stats':[res.statistic], 'step1_p_val':[res.pvalue]}))

    step1_res_table = pd.concat(calres_list, axis=0, ignore_index=True)
    step1_res_table["step1_p_adj"] = scipy.stats.false_discovery_control(step1_res_table.step1_p_val)
    step1_res_table = step1_res_table[step1_res_table["step1_p_adj"] < p_val]

    # ===== Step 2 =====
    # 與 whole edge 整合, 為了統合 edge ID
    # 知道 filtered network edge 在 background edge 的 ID
    tend_taxa2taxa_edge = background_edge[(background_edge['from'].isin(filter_edge['from'].to_list()))
                                        & (background_edge['to'].isin(filter_edge['to'].to_list()))]

    other_taxa2taxa_edge = background_edge[~background_edge['ID'].isin(tend_taxa2taxa_edge['ID'].to_list())]

    # browse each pathway & calculate
    calres_list2 = []
    total_pathway = list(path_koDict.keys())
    for current_path in total_pathway:
        current_taxa_list = path_taxaDict[current_path]

        # left col (in tendency network edge)
        # extract target edge
        tend_subedge = tend_taxa2taxa_edge[tend_taxa2taxa_edge["from"].isin(current_taxa_list) & tend_taxa2taxa_edge["to"].isin(current_taxa_list)]
        # store values
        l1, l2 = len(tend_subedge), len(tend_taxa2taxa_edge) - len(tend_subedge)

        # right col (the other edges in whole network)
        other_subedge = other_taxa2taxa_edge[other_taxa2taxa_edge["from"].isin(current_taxa_list) & other_taxa2taxa_edge["to"].isin(current_taxa_list)]
        r1, r2 = len(other_subedge), len(other_taxa2taxa_edge) - len(other_subedge)

        # fisher's exact test
        m = [[l1, r1],[l2, r2]]
        res = fisher_exact(m, alternative='greater')
        calres_list2.append(pd.DataFrame({'pathway':[current_path], 'step2_stats':[res.statistic], 'step2_p_val':[res.pvalue]}))

        step2_res_table = pd.concat(calres_list2, axis=0, ignore_index=True)
        step2_res_table["step2_p_adj"] = scipy.stats.false_discovery_control(step2_res_table.step2_p_val)
        step2_res_table  = step2_res_table [step2_res_table ["step2_p_adj"] < p_val]

    # ===== intergrate step1 and step2 =====
    both_sig_path = [path for path in step2_res_table["pathway"].to_list() if path in step1_res_table["pathway"].to_list()]
    # layout
    tmp1 = step1_res_table[step1_res_table["pathway"].isin(both_sig_path)]
    tmp2 = step2_res_table[step2_res_table["pathway"].isin(both_sig_path)]
    final_res_df = pd.concat([tmp1, tmp2.iloc[:,1:4]], axis=1)
    annotate = []
    for path in final_res_df["pathway"].values:
        map_name = pathwayInfo["pathway"][pathwayInfo["mapid"] == path].values[0]
        annotate.append(map_name)
    final_res_df["annotation"] = annotate
    return final_res_df


# ===== Main =====
background_edge = None
NC_dir = os.path.join(work_dir, "network_construction")
# check if background edge exists in tmp_files. Do not need to calculate again
if not os.path.isfile(os.path.join(NC_dir, "whole_taxa2taxa_edge.tsv")):
    make_whole_edge(num_processes=threads) # generate background edge table
    background_edge = pd.read_csv(os.path.join(NC_dir, "whole_taxa2taxa_edge.tsv"), sep="\t")
else:
    background_edge = pd.read_csv(os.path.join(NC_dir, "whole_taxa2taxa_edge.tsv"), sep="\t")

filter_dir = os.path.join(work_dir, "network_filter", target_label)
filter_edge = pd.read_csv(os.path.join(filter_dir, f"{target_label}_edge_f.csv"))
filter_node = pd.read_csv(os.path.join(filter_dir, f"{target_label}_node_f.csv"))
final_res = vital_pathway(filter_edge, filter_node, background_edge)
final_res.to_csv(os.path.join(tmp_dir, f"{target_label}_pathway.csv"), sep=",", index=False)
print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - pathway analysis done!")
sys.exit(0)
