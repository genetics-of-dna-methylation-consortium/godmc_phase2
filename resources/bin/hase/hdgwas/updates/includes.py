import time
import pandas as pd

class process_meth():
    """
    This makes reduced the time needed for reading in the methylation data
    """

    def __init__(self, IN_FILE):
        self.IN_FILE = IN_FILE

    def read_meth(self):
        if ((self.IN_FILE).split("/")[-1].split(".")[0] == "methylation_data"):
            start = time.time()
            chunk = pd.read_csv(self.IN_FILE, chunksize=100000, index_col=None)
            pd_df = pd.concat(chunk)

            # for i in ['\t', ' ']:
            #     chunk = pd.read_csv(self.IN_FILE, sep=i, chunksize=100000, index_col=None)
            #     pd_df = pd.concat(chunk)
            #     if pd_df.shape[1] > 1:
            #        break

            # if ((self.IN_FILE).split(".")[-1] == "csv"):    
            #     start = time.time()
            #     #Read the csv  data
            #     chunk = pd.read_csv(self.IN_FILE, chunksize=100000, index_col=None)
            #     pd_df = pd.concat(chunk)

            # if ((self.IN_FILE).split(".")[-1] == "txt"):    
            #     start = time.time()
            #     #Read the txt  data
            #     chunk = pd.read_csv(self.IN_FILE, sep='\t' chunksize=100000, index_col=None)
            #     pd_df = pd.concat(chunk) 
                
            #Transpose the DataFrame
            pd_df_transposed = pd_df.transpose().reset_index()
            pd_df_transposed.columns = pd_df_transposed.iloc[0]
            pd_df_transposed.drop(0, axis = 0, inplace = True) #initial index dropped 
            col_float = pd_df_transposed.iloc[:,1:].astype(float) #get float values
            iid = pd_df_transposed.iloc[:,0] #get object values - the columns
            df = pd.concat([iid, col_float], axis = 1) #merge to get the final df

            end = time.time()    
            print("Read data with chunks: ", (end-start),"sec")
            print(df.head())

        else:
            for i in ['\t', ' ']:
                df = pd.read_csv(self.IN_FILE, sep=i, index_col=None)
                if df.shape[1] > 1:
                   break
            else:
                raise ValueError(
                        'Cant read {} file; default settings: index=None, header=True; sep=tab or space '.format(file))
        
        return (df)