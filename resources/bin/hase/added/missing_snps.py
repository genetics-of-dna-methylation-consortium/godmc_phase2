import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Searches for missing SNPS after the SNP invertion step convertion")
parser.add_argument('-g', '--converted_folder', type = str, required=True, help='Path to the converted folder')
parser.add_argument('-n', '--study_name', type = str, required=True, help='Study name')

def search_missing_snps(DATA_PATH, STUDY_NAME):
    with pd.HDFStore(
        os.path.join(DATA_PATH, "probes", STUDY_NAME+"_original.h5"),'r') as f:
        data_original = f['probes']

    with pd.HDFStore(
        os.path.join(DATA_PATH, "probes", STUDY_NAME+".h5"),'r') as f:
        data = f['probes']

    snp_list = list(data["ID"])

    eq = data_original[~(data_original["ID"].isin(snp_list))]

    print(eq.head())
    
    #OUTPUT = "missing.csv"
    #eq.to_csv(os.path.join(DATA_PATH, OUTPUT))

if __name__ == "__main__":
    args = parser.parse_args()
    search_missing_snps(args.converted_folder, args.study_name)