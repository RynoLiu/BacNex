#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import time
import json
import os

def make_table_help():
    help = """
Usage: python make_table.py [options]

Parameters:
    -i,  --input_file    Path to the HUMAnN format table (ko.tsv).
    -o,  --outdir        The output directory.
    -h,  --help          Display the help message.
    """
    print(help)
    return

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate taxa-only relative abundance table based on HUMAnN ko.tsv')
parser.add_argument('-i', '--input', required=True, help='Path to the input HUMAnN ko.tsv file')
parser.add_argument('-o', '--outdir', required=True, help='The output directory')
args = parser.parse_args()
# parameters
input_file = args.input
outdir = args.outdir

# main function
def make_table(ko_table):
    """ The function generate the taxa-only relative abundance table based on HUMAnN ko.tsv
    Parameters
    ----------
        ko_table : HUMAnN ko.tsv

    Returns
    ----------
        taxa-only relative abundance table
        clean colnames ko.tsv

    """
    # clean id name from HUMAnN table
    o_id = ko_table.columns
    clean_id = [id.split("_merge_Abundance")[0] for id in o_id] # to rename final merged table

    # 1. taxa abundance table
    merge_taxa = {}
    for idx, c_val in ko_table[2:].iterrows(): # skip unmapped and ungrouped
        row_name = c_val.name.split('|') # list
        # skip ko id rows
        if len(row_name) <= 1:
            continue
        ko = row_name[0]
        taxa = row_name[1]
        val = c_val.values

        if taxa not in merge_taxa:
            # add new taxa keys
            merge_taxa[taxa] = val
        else:
            l_sum = np.sum([merge_taxa[taxa], val], axis=0).tolist() # conbine two values in list
            # update taxa values
            merge_taxa[taxa] = l_sum

    # make a new dataframe
    merged_table_taxa = pd.DataFrame(merge_taxa, index=clean_id)
    merged_table_taxa = merged_table_taxa.T
    merged_table_taxa = merged_table_taxa.loc[(merged_table_taxa!=0).any(axis=1)]
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - taxa table done!")

    # 2. clean colnames of ko.tsv
    rawdata = ko_table
    rawdata.rename(columns=dict(zip(rawdata.columns, clean_id)), inplace=True)
    rawdata = rawdata.T
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - functional table done!")

    return merged_table_taxa, rawdata


# taxa_koDict
def taxa_ko_dict(rawdata):
    # raw data : hum_table
    rawdata = rawdata.loc[:,rawdata.sum(0) > 0]
    row_val = rawdata.columns.to_list()

    # browse row
    taxa_koDict = {} # {taxa : [ko1, ko2, ..], ...}
    for w in row_val:
        if 'UNMAPPED' in w:
            continue
        if 'UNGROUPED' in w:
            continue
        if 'unclassified' in w:
            continue
        if(len(w.split('|')) == 1):
            continue
        ko, taxa = w.split('|')
        if taxa in taxa_koDict:
            taxa_koDict[taxa].append(ko)
        else:
            taxa_koDict[taxa] = [ko] # whole 的 taxa ko_list

    return taxa_koDict


# parser command
def run_make_table(options):
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - make_table running...")
    table = pd.read_csv(input_file, sep="\t", index_col=0) # ko.tsv
    taxa_data, hum_data = make_table(table)
    taxa_koDict = taxa_ko_dict(hum_data)
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - Saving files... (It may take a few minutes.)")
    taxa_data.to_csv(os.path.join(outdir, "taxa_table.csv"), sep=",", index=True)
    hum_data.to_csv(os.path.join(outdir, "func_table.tsv"), sep="\t", index=True)
    with open(os.path.join(outdir, 'taxa_koDict.json'), 'w') as fp1: json.dump(taxa_koDict, fp1)
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} - Done!")

if __name__ == "__main__":
    run_make_table(args)