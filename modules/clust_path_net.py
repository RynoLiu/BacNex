import pandas as pd
import numpy as np
import json
import sys
import os
import re

# parameters
work_dir = sys.argv[1] # taxa_table directory
target_pathway = sys.argv[2] # target pathway
target_label = sys.argv[3] # target cohort
cluster_id = sys.argv[4] # cluster id

# ===== prepare =====
script_dir = os.path.dirname(os.path.abspath(__file__))
tmp_dir = os.path.join(work_dir, "tmp_files")
# map_id to Kos
whole_komap_dir = os.path.join(script_dir, "whole_komap", "total_KO_pathway.txt")
pathwayInfo = pd.read_csv(whole_komap_dir, sep="\t")

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
path_taxaDict = {} # restore
for pk, pv in path_koDict.items():
    # browse taxa_koDict
    path_taxaDict[pk] = [tk for tk, tv in taxa_koDict.items() if set(pv).intersection(set(tv))]

# ===== main =====
cluster_method = re.split(r'(?=\d)', cluster_id)[0].rstrip("_")
edge_f = pd.read_csv(os.path.join(work_dir, "network_filter", target_label, f"{target_label}_edge_f.csv"))
node_f = pd.read_csv(os.path.join(work_dir, "network_filter", target_label, f"{target_label}_node_f.csv"))
clust_overview = pd.read_csv(os.path.join(work_dir, "network_analysis", f"{target_label}_cluster_overview", f"{target_label}_clust_node.csv"))
clust_subnet_node = clust_overview["id"][clust_overview[cluster_method] == cluster_id].to_list()

# Task 1 : get target pathway taxa and exists in target cluster
# output_1 : path_subnet_edge, path_subnet_node
edge_f["ID"] = range(edge_f.shape[0]) # add ID column
target_path_taxalist = path_taxaDict[target_pathway]
subnet_id = []
for i in range(len(edge_f)):
    row = edge_f.iloc[i,]
    if all( (row[key] in target_path_taxalist) and (row[key] in clust_subnet_node) for key in ["from", "to"] ):
        subnet_id.append(row.ID)
path_subnet_edge = edge_f[edge_f["ID"].isin(subnet_id)]
path_subnet_edge = path_subnet_edge.drop(columns=["ID"])
path_subnet_node = node_f[node_f["id"].isin(set(np.concatenate((path_subnet_edge["from"].values, path_subnet_edge["to"].values))))]

# Task 2 : target pahtway ko list from sub-network taxa
# ouput_2 : subnet_ko
target_mapko = path_koDict[target_pathway]
subnet_ko = set()
for i in path_subnet_node["id"].to_list():
    taxa_ko = [k for k in target_mapko if k in taxa_koDict[i]]
    for ko in taxa_ko:
        subnet_ko.add(ko)

# saved in tmp_files
path_subnet_edge.to_csv(os.path.join(tmp_dir, f"{target_label}_{cluster_id}_{target_pathway}_edge.csv"), index=False)
path_subnet_node.to_csv(os.path.join(tmp_dir, f"{target_label}_{cluster_id}_{target_pathway}_node.csv"), index=False)

# the KO set exist in taxa features and target pathway, used for KEGG pathway visualization
with open(os.path.join(tmp_dir, f"{target_label}_{cluster_id}_{target_pathway}_koset.txt"), "w") as fres:
    for ko in subnet_ko:
        fres.write(f"{ko}\n")