#!python
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import os

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description="This inputs the hash values from the converted probes to into the reference panel")
parser.add_argument('-i', '--input', type = str, required=True, help='Path to the initial reference panel')
parser.add_argument('-o', '--output', type = str, required=True, help='Path to save the hashed reference panel')
parser.add_argument('-t', '--hash-table', type = str, required=True, help='Path to the hash table of the converted probes')

def rf_hash(ref_panel, output_folder, hash_table):
    hash_df = pd.read_csv(hash_table, compression='gzip', sep ="\t")

    print("The hash table values...")
    print(hash_df.head)

    T_new = hash_df[hash_df["allele"] == "T"]["keys"].values[0]
    A_new = hash_df[hash_df["allele"] == "A"]["keys"].values[0]
    C_new = hash_df[hash_df["allele"] == "C"]["keys"].values[0]
    G_new = hash_df[hash_df["allele"] == "G"]["keys"].values[0]
    
    ref_original = pd.read_csv(ref_panel, compression='gzip', sep =" ")
    
    chunk = pd.read_csv(ref_panel, chunksize=100000000, compression='gzip', sep =" ")
    ref_original = pd.concat(chunk)

    print("The initial reference panel")
    print(ref_original.head)

    T_old = ref_original[ref_original["str_allele1"] == "T"].head(1)["allele1"].values[0]
    A_old = ref_original[ref_original["str_allele1"] == "A"].head(1)["allele1"].values[0]
    C_old = ref_original[ref_original["str_allele1"] == "C"].head(1)["allele1"].values[0]
    G_old = ref_original[ref_original["str_allele1"] == "G"].head(1)["allele1"].values[0]
    
    ref_original[["allele_1_computed", "allele_2_computed"]]= ""
    for i in range(len(ref_original)):
        if (ref_original.loc[i, "str_allele1"] == "T") & (ref_original.loc[i, "allele1"] == T_old):
            ref_original.loc[i, "allele_1_computed"] = T_new

        elif (ref_original.loc[i, "str_allele1"] == "A") & (ref_original.loc[i, "allele1"] == A_old):
            ref_original.loc[i, "allele_1_computed"] = A_new

        elif (ref_original.loc[i, "str_allele1"] == "C") & (ref_original.loc[i, "allele1"] == C_old):
            ref_original.loc[i, "allele_1_computed"] = C_new

        elif (ref_original.loc[i, "str_allele1"] == "G") & (ref_original.loc[i, "allele1"] == G_old):
            ref_original.loc[i, "allele_1_computed"] = G_new
        else:
             ref_original.loc[i, "allele_1_computed"] = ""

        #recomputation for allele 2 in the refence panel
        if (ref_original.loc[i, "str_allele2"] == "T") & (ref_original.loc[i, "allele2"] == T_old):
            ref_original.loc[i, "allele_2_computed"] = T_new

        elif (ref_original.loc[i, "str_allele2"] == "A") & (ref_original.loc[i, "allele2"] == A_old):
            ref_original.loc[i, "allele_2_computed"] = A_new

        elif (ref_original.loc[i, "str_allele2"] == "C") & (ref_original.loc[i, "allele2"] == C_old):
            ref_original.loc[i, "allele_2_computed"] = C_new

        elif (ref_original.loc[i, "str_allele2"] == "G") & (ref_original.loc[i, "allele2"] == G_old):
            ref_original.loc[i, "allele_2_computed"] = G_new
        else:
             ref_original.loc[i, "allele_2_computed"] = ""
    
    ref_original.drop(["allele1", "allele2"], axis = 1, inplace=True)
    ref_original.rename({"allele_1_computed":"allele1", "allele_2_computed":"allele2"}, axis = 1, inplace=True)
    
    ref_new = ref_original[["ID", "bp", "str_allele1", "str_allele2", "CHR", "allele1", "allele2"]]

    print("Hash values are now inputed into the reference panel")
    print(ref_new.head)

    ref_new.to_csv(os.path.join(output_folder, 'ref-hrc.ref.gz'), compression='gzip', sep =" ")
    print("Hashed reference panel in now saved")

if __name__ == "__main__":
    args = parser.parse_args()
    rf_hash(args.input, args.output, args.hash_table)

