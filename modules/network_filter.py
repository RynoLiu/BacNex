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
opt_PL_w = float(sys.argv[3]) # PL weight threshold
opt_edge_w = float(sys.argv[4]) # edge weight threshold
select_task = sys.argv[5]


# ===== prepare =====
script_dir = os.path.dirname(os.path.abspath(__file__))

# PreLect result
PLres = pd.read_csv(os.path.join(work_dir, "prelect_dir", "PLres.csv"), index_col=0)


# ===== main =====
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
edge = pd.read_csv(os.path.join(work_dir, "network_construction", target_label, f"{target_label}_edge.csv"))
node = pd.read_csv(os.path.join(work_dir, "network_construction", target_label, f"{target_label}_node.csv"))
# filtered by weight
edge_f = edge[
(edge["from"].isin(PL_fl)) &
(edge["to"].isin(PL_fl)) &
(edge["weight"] > opt_edge_w)]

node_f = node[node["id"].isin(PL_fl)]


# save
edge_f.to_csv(os.path.join(work_dir, "tmp_files", f"{target_label}_edge_f.csv"), index=False)
node_f.to_csv(os.path.join(work_dir, "tmp_files", f"{target_label}_node_f.csv"), index=False)
print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - Network filter done!")
