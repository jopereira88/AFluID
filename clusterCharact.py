#!/usr/bin/python3
''' 
Part of the pipeline installation: creates clustes annotation files.
Expects unchanged filenames and dirtree, if they are changed they will 
need to be changed in the scripts
'''

from flu_utils import metadata_dict,parse_clstr
from collections import Counter
import pickle

meta=metadata_dict('metadata/flu_metadata_4.1_flt.csv',[9,10,12])
clusters=parse_clstr('cluster_db/infDNAClusters.clstr')

cluster_table={}
for key in clusters:
    cluster_table[key]=[Counter(),Counter(),Counter()]
    for value in clusters[key]:
        cluster_table[key][0][meta[value][0]]+=1
        cluster_table[key][1][meta[value][1]]+=1
        cluster_table[key][2][meta[value][2]]+=1
for elem in cluster_table:
    cluster_table[elem][0]=dict(cluster_table[elem][0])
    cluster_table[elem][1]=dict(cluster_table[elem][1])
    cluster_table[elem][2]=dict(cluster_table[elem][2])

clusters=parse_clstr('cluster_db/infDNAClusters.clstr',access_only=False)
cluster_reps={}
for cl in clusters:
    rep=[item for item in clusters[cl] if '*' in item]
    cluster_reps[cl]=rep[0].split(' ')[1].replace('...','').replace('>','')
for cl in cluster_table:
    cluster_table[cl].append(cluster_reps[cl])

with open('cluster_db/dict_clusters.pkl','wb') as save:
    pickle.dump(cluster_table,save)
with open('cluster_db/cluster_desc.txt','w') as file:
    file.write(f'Cluster\tRepresentative\tGenotypes\tSegments\tHosts\n')
    for elem in cluster_table:
        file.write(f'{elem}\t{cluster_table[elem][3]}\t{cluster_table[elem][0]}\t{cluster_table[elem][1]}\t{cluster_table[elem][2]}\n')