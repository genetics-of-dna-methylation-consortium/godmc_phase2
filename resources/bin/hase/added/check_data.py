import os
import argparse
import pandas as pd
import h5py

parser = argparse.ArgumentParser(description="Checks the data after convertion")
parser.add_argument('-g', '--converted_folder', type = str, required=True, help='Path to the converted folder')
parser.add_argument('-n', '--study_name', type = str, required=True, help='Study name')

def check_data(DATA_PATH, STUDY_NAME):
    l = os.listdir(os.path.join(DATA_PATH, 'genotype'))
    l.sort()
    l.sort(key = len)
    f = 0

    for i in range(len(l)):
        a, b = h5py.File(os.path.join(DATA_PATH, 'genotype', l[i]), 'r')['genotype'][...].shape
        f += a

    print("Size of Genotype Data: ", f)

    with pd.HDFStore(os.path.join(DATA_PATH, "probes", STUDY_NAME+".h5"),'r') as f:
        df = f['probes']
        
    print("Size of Probes: ", len(df))

if __name__ == "__main__":
    args = parser.parse_args()
    check_data(args.converted_folder, args.study_name)