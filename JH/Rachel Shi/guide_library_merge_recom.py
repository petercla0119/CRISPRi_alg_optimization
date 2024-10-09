#Duplicates guide library rows to account for recombinant reads, removes G in keys, adds "R_" to beginning of sgRNA name for recom rows
import pandas as pd
import csv
def recom_first_read(seq):
    return seq[0:19]*2
def recom_name_change(name):
    return 'R_' + name
df = pd.read_csv('20200513_library_1_2_unbalanced_dJR051.csv')
#df = pd.read_csv('guide_library_test.csv')
df['protospacer_A'] = df['protospacer_A'].str.slice(start=1,stop=20)
df['protospacer_B'] = df['protospacer_B'].str.slice(start=1,stop=20)
df['protospacer_A'] = df['protospacer_A'] + df['protospacer_B']
recom = df.copy()
recom['protospacer_A'] = recom['protospacer_A'].apply(recom_first_read)
recom['protospacer_A'] = recom['protospacer_A'].apply(recom_first_read)
recom['sgID_A'] = recom['sgID_A'].apply(recom_name_change)
recom['sgID_B'] = recom['sgID_B'].apply(recom_name_change)
frames = [df,recom]
result = pd.concat(frames, ignore_index=True)
result.to_csv('guide_library_merging_recom.csv')