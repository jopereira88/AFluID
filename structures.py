segment_syns={'1':['PB2','segment 1', 'segment: 1','chromossome: 1'],'2':['PB1','segment 2','PB1-F2','segment: 2'],'3':['PA','segment 3','PA-X','segment: 3'],
              '4':['hemagglutinin','haemagglutinin','HA1','HA2','HA','segment 4','segment: 4'],'5':['NP','segment 5','nucleocapside protein','segment: 5'],
              '6':['neuraminidase', 'NA','segment 6','segment: 6'],'7':['segment 7','M1','M2','matrix protein','segment: 7'],
              '8':['segment 8', 'NS1' ,'NS2', 'NEP', 'nonstructural protein', 'nuclear export protein','segment: 8']}
iupac_to_int={'PB2':1,'PB1':2,'PA':3,'HA':4,'NP':5,'NA':6,'MP':7,'NS':8}
flagdict={"Fasta Preprocess":{"Rejected Sequences":[]},"Master":{"C_BLAST":True,"L-BLAST":True,"flumut":True,"genin":True,"nextclade":True,"getref":False,"single":True}, 
          "CD-HIT":{"Unclustered Sequences":[]},
          "BLAST":{"Sequences unassigned against cluster representatives":[],"Sequences unassigned against local database":[]},
          'Sample':{'HA':[],'NA':[],'PB2':[],'PB1':[],'PA':[],'NP':[],'MP':[],'NS':[],'Genotype':'','HA_ref':[],'NA_ref':[],
                    'PB2_ref':[],'PB1_ref':[], 'PA_ref':[], 'NP_ref':[],'MP_ref':[],'NS_ref':[],'PB2_len':0,'PB1_len':0,
                    'PA_len':0,'HA_len':0,'NP_len':0,'MP_len':0,'NA_len':0,'NS_len':0,'H_gen':'','N_gen':''},
          "Final Report":{"Sequences for FluMut":[],"Sequences for NextClade":{'H1':[],'H3':[],'H5':[]},"Sequences for GenIn":[],'Get References':[]}}
taxa_dict={}
geo_dict={}
