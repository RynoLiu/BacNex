import pandas as pd
import numpy as np
import json
import time
import copy
import sys
import os


# parameters
work_dir = sys.argv[1] # taxa_table directory
target_label = str(sys.argv[2]) # target label
size_opt = sys.argv[3] # size options
opt_PL_w = float(sys.argv[4]) # PL weight threshold
opt_edge_w = float(sys.argv[5]) # edge weight threshold
select_rank = sys.argv[6] # select rank

# visualization refined options
size_range_min = float(sys.argv[7]) # size range min
size_range_max = float(sys.argv[8]) # size range max
width_range_min = float(sys.argv[9]) # width range min
width_range_max = float(sys.argv[10]) # width range max

select_task = sys.argv[11] # selected classification task


# ===== prepare =====
script_dir = os.path.dirname(os.path.abspath(__file__))
# rank = ["phylum", "class", "order", "family", "genus", "none"]
sp_map = pd.read_csv(os.path.join(script_dir, "taxa_lineage", "species_map.csv")) # group info map

# PreLect result
PLres = pd.read_csv(os.path.join(work_dir, "prelect_dir", "PLres.csv"), index_col=0)

# MC FC nessesary table
taxa_table = pd.read_csv(os.path.join(work_dir, "taxa_table.csv"), index_col=0)
meta = pd.read_csv(os.path.join(work_dir, "tmp_files", "meta.csv"), index_col=0)


# ===== Main functions =====
# node table
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
            group = sp_map[select_rank][sp_map["genus"] == genus_name].values[0]
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


# ===== main =====
# taxa_koDict
taxa_koDict = None
with open(os.path.join(work_dir, "taxa_koDict.json"), 'r') as fp: taxa_koDict = json.load(fp)

# PLlist
PLres = pd.read_csv(os.path.join(work_dir, "prelect_dir", "PLres.csv"), index_col=0)
# extract tendency taxa
tend_PL = None
PL_fl = None

if select_task in ["BC", "Rg", "Cox"]:
    if target_label != "All":
        tend_PL = PLres[PLres["tendency"] == target_label].iloc[:, 0:3]
        PL_fl = tend_PL["FeatName"][abs(tend_PL["coef"]) >= opt_PL_w].to_list() # feature list
    else:
        tend_PL = PLres[PLres["selected"] == "Selected"].iloc[:, 0:3]
        PL_fl = tend_PL["FeatName"][abs(tend_PL["coef"]) >= opt_PL_w].to_list() # feature list

elif select_task == "MC":
    tend_col = f"tendency_{target_label}"
    coef_col = f"coef_{target_label}"
    tend_PL = PLres[PLres[tend_col] == target_label].loc[:, ["FeatName", coef_col, tend_col]]
    PL_fl = tend_PL["FeatName"][abs(tend_PL[coef_col]) >= opt_PL_w].to_list() # feature list

# edge weight filter
edge = pd.read_csv(os.path.join(work_dir, "network_construction", f"{target_label}_edge.csv"))
# filtered by weight
edge_f = edge[
(edge["from"].isin(PL_fl)) &
(edge["to"].isin(PL_fl)) &
(edge["weight"] > opt_edge_w)]

# node table of edge_f
node_f = node_table(edge_f, size_opt)
# scale width of edge_f
edge_f = width_scaler(edge_f)


# save
edge_f.to_csv(os.path.join(work_dir, "tmp_files", f"{target_label}_edge_f.csv"), index=False)
node_f.to_csv(os.path.join(work_dir, "tmp_files", f"{target_label}_node_f.csv"), index=False)
print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - Network filter done!")
