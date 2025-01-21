#!python
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import os
import sys


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

    # To avoid behaviour different to original loop, check for any instances where only one allele is NA in both dataframes
    # These would also have been missed by the original loop resulting in errors when type casting allele1/2 below
    any_x_na_mismatches = (merged["allele1_x"].isna() ^ merged["allele2_x"].isna()).any()
    any_y_na_mismatches = (merged["allele1_y"].isna() ^ merged["allele2_y"].isna()).any()

    # Error and report if the reference or input file have probes with one missing ref/alt allele
    if any_x_na_mismatches:
        input_reference_path = os.path.join("./resources/bin/hase/data", "ref-hrc.ref.gz")
        sys.exit("Error: Partially missing alleles in reference panel file. Please check the input file: {}".format(input_reference_path))
    elif any_y_na_mismatches:
        sys.exit("Error: Partially missing alleles in probes input file. Please check the input file: {}".format(new_original))

    # Provided there are no partially missing ref/alt alleles, use combine first to fill NAs for alleles not in the reference panel
    # This caters to SNPS that are not in the reference panel
    merged["allele1"] = merged["allele1_x"].combine_first(merged["allele1_y"])
    merged["allele2"] = merged["allele2_x"].combine_first(merged["allele2_y"])

    merged.rename({"CHR_y":"CHR", "bp_y":"bp"}, axis = 1, inplace=True)
    merged = merged.reindex(columns=["CHR", "ID", "distance", "bp", "allele1", "allele2"])

    merged["ID"] = merged["ID"].apply(lambda x: (str(x)))
    merged["allele1"] =  merged["allele1"].apply(lambda x: (int(x)))
    merged["allele2"] =  merged["allele2"].apply(lambda x: (int(x)))

    merged.to_hdf(os.path.join(probes_folder, project_name+".h5"), key='probes',
                format='table',data_columns=True, append=True, complib='zlib',complevel=9, min_itemsize = 45)

    print("Script ran sucessfully")

if __name__ == "__main__":
    args = parser.parse_args()
    inv_probes(args.probes_folder, args.project_name)
