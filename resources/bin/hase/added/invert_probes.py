#!python
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import os


pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description="This inputs inverts the alleles in case the columns are wrong")
parser.add_argument('-f', '--probes_folder', type = str, required=True, help='Path to the probes folder')
parser.add_argument('-n', '--project_name', type = str, required=True, help='Name of the project')

def inv_probes(probes_folder, project_name):
    old_original = os.path.join(probes_folder, project_name + ".h5")
    new_original = os.path.join(probes_folder, project_name + "_original.h5")
    os.rename(old_original, new_original)
    
    chunk = pd.read_hdf(os.path.join(probes_folder, project_name + "_original.h5"), 'probes', chunksize=100000, compression='gzip', sep =" ")
    df = pd.concat(chunk)
    
    ref_HRC_chunked = pd.read_csv(os.path.join("./resources/bin/hase/data", "ref-hrc.ref.gz"), compression='gzip', sep =" ", chunksize=100000)
    ref_HRC = pd.concat(ref_HRC_chunked)

    merged = pd.merge(ref_HRC, df, on="ID", how = "right")
    flipped_alleles = merged[merged["allele1_x"] != merged["allele1_y"]]
    print ("{} Flipped Alleles Detected".format((len(flipped_alleles))))

    merged["allele1"] = merged["allele1_x"]
    merged["allele2"] = merged["allele2_x"]

    merged.rename({"CHR_y":"CHR", "bp_y":"bp"}, axis = 1, inplace=True)
    merged = merged.reindex(columns=["CHR", "ID", "distance", "bp", "allele1", "allele2"])

    #This caters dor SNPS that are not in the reference panel
    for i in range(len(merged)):
        if (pd.isnull(merged.loc[i, "allele1"])) & (pd.isnull(merged.loc[i, "allele2"])):
            x = df[df["ID"] ==  merged.loc[i, "ID"]]["allele1"].values[0]
            y =df[df["ID"] ==  merged.loc[i, "ID"]]["allele2"].values[0]
            merged.loc[i, "allele1"] =  x
            merged.loc[i, "allele2"] =  y

    merged["ID"] = merged["ID"].apply(lambda x: (str(x)))
    merged["allele1"] =  merged["allele1"].apply(lambda x: (int(x)))
    merged["allele2"] =  merged["allele2"].apply(lambda x: (int(x)))

    merged.to_hdf(os.path.join(probes_folder, project_name+".h5"), key='probes',
                format='table',data_columns=True, append=True, complib='zlib',complevel=9, min_itemsize = 45)

    print("Script ran sucessfully")

if __name__ == "__main__":
    args = parser.parse_args()
    inv_probes(args.probes_folder, args.project_name)