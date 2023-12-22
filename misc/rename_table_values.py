import sys
import pandas as pd

fname = sys.argv[1]
cols = sys.argv[2]
conversion_fname = sys.argv[3]
key_col = int(sys.argv[4])
value_col = int(sys.argv[5])

data = pd.read_csv(fname, sep = "\t", header = None) 
sample_conversion = pd.read_csv(conversion_fname, sep = "\t") 
dict_conversion = dict(zip(sample_conversion.iloc[:,key_col], sample_conversion.iloc[:,value_col]))

for col in map(int,cols.split(',')):
    id_list = list(data.iloc[:,col])
    data.iloc[:,col] = [dict_conversion.get(id, id) for id in id_list]


data.to_csv(sys.stdout, sep="\t", na_rep="NA", header = False, index = False)