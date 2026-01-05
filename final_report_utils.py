from collections import defaultdict
import ast
from structures import geo_dict, taxa_dict, muts_interest, seg_lens, muts_loci_meaning
import pandas as pd
import os

def rollup_geo(location, geo_dict, mode='region'):
    if location in geo_dict.keys():
        if mode=='region':
            return geo_dict[location][0]
        elif mode=='continent':
            return geo_dict[location][1]
    else:
        return 'Unknown'

def search_tax_level(term, taxa_dict):
    tax_level={0:'Species', 1:'Genus', 2:'Subfamily', 3:'Family', 4:'Division', 5:'Order', 6:'Class'}
    for key in taxa_dict:
        for level in range(len(taxa_dict[key])):
            if term == taxa_dict[key][level]:
                return tax_level[level] , key
    return 'Unknown'

def rollup_taxa(term, taxa_dict, level='Genus'):
    level_tax={'SPECIES':0, 'GENUS':1, 'SUBFAMILY':2, 'FAMILY':3, 'DIVISION':4, 'ORDER':5, 'CLASS':6}
    t_lev,key=search_tax_level(term, taxa_dict)
    if t_lev.upper() == level.upper():
        return taxa_dict[key][level_tax[level.upper()]]
    elif level_tax[t_lev.upper()] < level_tax[level.upper()]:
        return taxa_dict[key][level_tax[level.upper()]]
    elif level_tax[t_lev.upper()] > level_tax[level.upper()]:
        return taxa_dict[key][level_tax[t_lev.upper()]]
    else:
        return 'Unknown'
    
def _rollup_counts(counter, mapper, drop_unknown=False):
    """
    Generic aggregator: maps each key via `mapper(key) -> new_key` and sums values.
    If `drop_unknown` is True, entries that map to 'Unknown' are discarded.
    """
    out = defaultdict(float)
    for k, v in counter.items():
        new_k = mapper(k)
        if new_k == 'Unknown' and drop_unknown:
            continue
        out[new_k] += v
    # cast back to int if all inputs were ints and sums are integral
    if all(isinstance(v, int) for v in counter.values()):
        return {k: int(v) for k, v in out.items()}
    return dict(out)

def rollup_counts_geo(counter, geo_dict, mode='region', drop_unknown=False):
    """
    Roll up location counts to ISO region or continent using your `rollup_geo`.
      - mode: 'region' or 'continent'
    """
    return _rollup_counts(counter, lambda loc: rollup_geo(loc, geo_dict, mode=mode),
                          drop_unknown=drop_unknown)

def rollup_counts_taxa(counter, taxa_dict, level='Genus', drop_unknown=False):
    """
    Roll up taxonomic terms to the target `level` using your `rollup_taxa`.
      - level options (case-insensitive): Species, Genus, Subfamily, Family, Division, Order, Class
    """
    return _rollup_counts(counter, lambda term: rollup_taxa(term, taxa_dict, level=level),
                          drop_unknown=drop_unknown)

def normalize(counter):
    """
    Convert absolute counts to proportions that sum to 1.0.
    """
    total = sum(counter.values())
    if total == 0:
        return {k: 0.0 for k in counter}
    return {k: v / total for k, v in counter.items()}

def clean_host(host_cell, taxa_dict, level='Species', drop_unknown=True):
    host_dict=host_cell.strip()
    host_dict=host_dict.replace('%','')
    host_dict=ast.literal_eval(host_dict)
    for key in host_dict:
        host_dict[key]=float(host_dict[key])
    host_dict=rollup_counts_taxa(host_dict, taxa_dict, level, drop_unknown)
    if round(sum(host_dict.values()),3)>1.1:
        for key in host_dict:
            host_dict[key]=round(host_dict[key]/100,3)
    else:
        for key in host_dict:
            host_dict[key]=round(host_dict[key],3)
    return host_dict

def clean_country(country_cell, geo_dict, level='country', drop_unknown=True):
    country_dict=country_cell.strip()
    country_dict=country_dict.replace('%','')
    country_dict=ast.literal_eval(country_dict)
    for key in country_dict:
        country_dict[key]=float(country_dict[key])
    country_dict=rollup_counts_geo(country_dict, geo_dict, level, drop_unknown)
    if sum(country_dict.values())>1.0:
        for key in country_dict:
            country_dict[key]=round(country_dict[key]/100,3)
    else:
        for key in country_dict:
            country_dict[key]=round(country_dict[key],3)
    return country_dict

def clean_genotype(genotype_cell):
    genotype_dict=genotype_cell.strip()
    genotype_dict=genotype_dict.replace('%','')
    genotype_dict=ast.literal_eval(genotype_dict)
    for key in genotype_dict:
        genotype_dict[key]=float(genotype_dict[key])
    if sum(genotype_dict.values())>1.0:
        for key in genotype_dict:
            genotype_dict[key]=round(genotype_dict[key]/100,3)
    else:
        for key in genotype_dict:
            genotype_dict[key]=round(genotype_dict[key],3)
    return genotype_dict

def convert_muts(flags,muts_loci_meaning):
    """
    Convert mutations from flags to a list of dictionaries with detailed info.
    """
    muts={}
    for segment in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
        if flags['Sample'][f'{segment}_muts']!=[]:
            for mut in flags['Sample'][f'{segment}_muts']:
                muts[muts_loci_meaning[mut][0]]=muts_loci_meaning[mut][1]

    return muts

def create_sample_card(flags):
    """
    Create HTML card for Sample tab in final report.
    """
    GENOTYPE=flags['Sample']['Genotype']
    H=flags['Sample']['H_gen']
    N=flags['Sample']['N_gen']
    CLADE=flags['Sample']['clade'] if flags['Sample']['clade'] else 'N/A'
    C_BLAST=flags['Master']['C_BLAST']
    L_BLAST=flags['Master']['L_BLAST']
    NEXTCLADE=flags['Master']['nextclade']
    FLUMUT=flags['Master']['flumut']
    GENIN=flags['Master']['genin']
    SINGLE=flags['Master']['single']
    card=f"""
    <div class="card" id="sample-card">
  <p>
    <strong>Genotype:</strong> {GENOTYPE}
    &nbsp;·&nbsp; <strong>H:</strong> {H}
    &nbsp;·&nbsp; <strong>N:</strong> {N}
    &nbsp;·&nbsp; <strong>Clade:</strong> {CLADE}
  </p>

  <div class="chips" id="pipeline-chips">
    <!-- Preencha com chips de estado dos módulos da pipeline -->
    <!-- Exemplo (repita/ajuste conforme preciso): -->
    <span class="chip">C_BLAST: <b>{C_BLAST}</b></span>
    <span class="chip">L_BLAST: <b>{L_BLAST}</b></span>
    <span class="chip">Nextclade: <b>{NEXTCLADE}</b></span>
    <span class="chip">FluMut: <b>{FLUMUT}</b></span>
    <span class="chip">GenIn: <b>{GENIN}</b></span>
    <span class="chip">Single: <b>{SINGLE}</b></span>
  </div>
</div>"""
    return card

def create_sample_table(df_segments,flags,seg_lens,mappings):
    """
    Create HTML table for Sample tab in final report.
    """
    #seg lengths thresholds based on seg_lens from structures.py
    len_thres={
        'PB2': (seg_lens['PB2_min'], seg_lens['PB2_max']),
        'PB1': (seg_lens['PB1_min'], seg_lens['PB1_max']),
        'PA': (seg_lens['PA_min'], seg_lens['PA_max']),
        'HA': (seg_lens['HA_min'], seg_lens['HA_max']),
        'NP': (seg_lens['NP_min'], seg_lens['NP_max']),
        'NA': (seg_lens['NA_min'], seg_lens['NA_max']),
        'MP': (seg_lens['MP_min'], seg_lens['MP_max']),
        'NS': (seg_lens['NS_min'], seg_lens['NS_max'])
    }
    reverse_map={v: k for k, v in mappings.items()}
    seq_names={}
    for segment in len_thres.keys():
        for name in flags['Sample'][segment]:
            seq_names[mappings[name]]=segment
    table='''
    <table class="tabela" id="sample-segments">
    <thead>
    <tr>
        <th>Segment</th>
        <th>Length</th>
        <th>Reference</th>
        <th>Mutations</th>
        <th>Note</th>
    </tr>
    </thead>
    <tbody>
    '''
    # Create table lines
    df=pd.read_table(df_segments,index_col=False)
    for row in range(len(df)):
        seq=df.loc[row,'SAMPLE_NAME']
        segment=seq_names[seq] if seq in seq_names else "Not Available"
        length=flags['Sample'][f'{segment}_len']
        if segment in len_thres:
            min_len, max_len = len_thres[segment]
            if length < min_len:
                length_class = "len-below"
                note='Length outside of expected range'
            elif length > max_len:
                length_class = "len-above"
                note='Length outside of expected range'
            else:
                length_class = "len-ok"
                note=''
            length_html = f'<span class="{length_class}">{length}</span>'
        else:
            length_html = str(length)
        reference=flags['Sample'][f'{segment}_ref'][0] if flags['Sample'][f'{segment}_ref'] else "Not Available"
        mutations=''.join(flags['Sample'][f'{segment}_muts']) if flags['Sample'][f'{segment}_muts'] else "None"
        if note!="":
            note+="; Segment not found" if segment=="Not Available" else ""
        else:
            note="Segment not found" if segment=="Not Available" else ""
        table += f"""
            <tr>
              <td>{segment}</td>
              <td>{length_html}</td>
              <td>{reference}</td>
              <td>{mutations}</td>
              <td>{note}</td>
            </tr>
        """
    # Close the table
    table += """
      </tbody>
    </table>
    """
    return table
def create_segment_background_table(df_segments,flags,seg_lens,mappings):
    """
    Placeholder function for creating Segment Background table.
    """
    #seg lengths thresholds based on seg_lens from structures.py
    len_thres={
        'PB2': (seg_lens['PB2_min'], seg_lens['PB2_max']),
        'PB1': (seg_lens['PB1_min'], seg_lens['PB1_max']),
        'PA': (seg_lens['PA_min'], seg_lens['PA_max']),
        'HA': (seg_lens['HA_min'], seg_lens['HA_max']),
        'NP': (seg_lens['NP_min'], seg_lens['NP_max']),
        'NA': (seg_lens['NA_min'], seg_lens['NA_max']),
        'MP': (seg_lens['MP_min'], seg_lens['MP_max']),
        'NS': (seg_lens['NS_min'], seg_lens['NS_max'])
    }
    reverse_map={v: k for k, v in mappings.items()}
    seq_names={}
    for segment in len_thres.keys():
        for name in flags['Sample'][segment]:
            seq_names[mappings[name]]=segment
    table = """
    <div class="sb-table-wrap">
    <table class="tabela" id="sb-flat-table">
    <thead>
      <tr>
        <th>Segment</th>
        <th>Length</th>
        <th>Reference</th>
        <th>Host</th>
        <th>Geography</th>
        <th>%ID</th>
        <th>Cluster</th>
        <th>Genotype</th>
        <th>Assigned_By</th>
        <th>Note</th>
      </tr>
    </thead>
    <tbody>"""
    df=pd.read_table(df_segments,index_col=False)
    mask=df['ASSIGNED_BY']!='L-BLAST'
    df.loc[mask,'HOST']=df.loc[mask,'HOST'].apply(lambda x: clean_host(x,taxa_dict,level='Species',drop_unknown=True))
    df.loc[mask,'GENOTYPE']=df.loc[mask,'GENOTYPE'].apply(lambda x: clean_genotype(x))
    df.loc[mask,'COUNTRY']=df.loc[mask,'COUNTRY'].apply(lambda x: ast.literal_eval(x) if x.startswith('{') else x)
    for row in range(len(df)):
        seq=df.loc[row,'SAMPLE_NAME']
        segment=seq_names[seq] if seq in seq_names else "Not Available"
        length=flags['Sample'][f'{segment}_len']
        if segment in len_thres:
            min_len, max_len = len_thres[segment]
            if length < min_len:
                length_class = "len-below"
                note='Length outside of expected range'
            elif length > max_len:
                length_class = "len-above"
                note='Length outside of expected range'
            else:
                length_class = "len-ok"
                note=''
            length_html = f'<span class="{length_class}">{length}</span>'
        else:
            length_html = str(length)
        reference=flags['Sample'][f'{segment}_ref'][0] if flags['Sample'][f'{segment}_ref'] else "Not Available"
        host_species=df.loc[row,'HOST']
        geo_countries=df.loc[row,'COUNTRY']
        if mask.loc[row]:
          host_order=rollup_counts_taxa(host_species,taxa_dict,level='Order',drop_unknown=True)
          host_class=rollup_counts_taxa(host_species,taxa_dict,level='Class',drop_unknown=True)
          geo_regions=rollup_counts_geo(geo_countries,geo_dict,mode='region',drop_unknown=True)
          geo_continents=rollup_counts_geo(geo_countries,geo_dict,mode='continent',drop_unknown=True)
        else:
          host_order=rollup_taxa(host_species,taxa_dict,level='Order')
          host_class=rollup_taxa(host_species,taxa_dict,level='Class')  
          geo_regions=rollup_geo(geo_countries,geo_dict,mode='region')
          geo_continents=rollup_geo(geo_countries,geo_dict,mode='continent')
        perc_id=df.loc[row,'%ID']
        cluster=df.loc[row,'CLUSTER']
        genotype=df.loc[row,'GENOTYPE']
        assigned_by=df.loc[row,'ASSIGNED_BY']
        if note!="":
            note+="; Segment not found" if segment=="Not Available" else ""
        else:
            note="Segment not found" if segment=="Not Available" else ""
        table += f"""
        <tr>
        <td>{segment}</td>
        <td>{length_html}</td>
        <td>{reference}</td>
        <td>
          <span class="host host--species">{host_species}</span>
          <span class="host host--order">{host_order}</span>
          <span class="host host--class">{host_class}</span>
        </td>
        <td>
          <span class="geo geo--country">{geo_countries}</span>
          <span class="geo geo--region">{geo_regions}</span>
          <span class="geo geo--continent">{geo_continents}</span>
        </td>
        <td>{perc_id}</td>
        <td>{cluster}</td>
        <td>{genotype}</td>
        <td>{assigned_by}</td>
        <td>{note}</td>
      </tr>  
        """
    table += """    
        </tbody>
        </table>
        </div>
            """
    return table
def create_segment_mutations_table(flags,muts_loci_meaning):
    """
    Placeholder function for creating Segment Mutations table.
    """
    muts=convert_muts(flags,muts_loci_meaning)
    table = """
    <table class="tabela" id="mutations-table">
    <thead>
      <tr>
        <th>Mutation</th>
        <th>Biological effects</th>
      </tr>
    </thead>
    <tbody>"""
    for mut in muts:
        meaning=muts[mut]
        table += f"""
        <tr>
        <td>{mut}</td>
        <td>{meaning}</td>
      </tr>  
        """
    table += """    
        </tbody>
        </table>
            """
    return table
    
html_skeleton = r"""<!DOCTYPE html>
<html lang="EN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title><!--FILENAME--> Final Report </title>

  <style>
    /* =========================
       Root tokens & base
       ========================= */
    :root{
      --bg:#fff; --fg:#111827; --muted:#6b7280;
      --accent:#2563eb; --accent-weak:#e5edff;
      --border:#e5e7eb; --card:#f9fafb;
      --sans: system-ui,-apple-system, Segoe UI, Roboto, Noto Sans, Ubuntu, Cantarell, Arial;
    }
    html,body{height:100%}
    body{margin:0;font-family:var(--sans);color:var(--fg);background:var(--bg);line-height:1.55}

    header{padding:1rem;border-bottom:1px solid var(--border)}
    .container{max-width:1100px;margin:0 auto;padding:1rem}

    /* =========================
       Tabs (no-JS)
       ========================= */
    .tabs{border:1px solid var(--border); border-radius:.75rem; background:#fff; overflow:hidden}
    .tab-controls{display:flex; gap:.5rem; padding:.5rem; border-bottom:1px solid var(--border); background:var(--card); flex-wrap:wrap}
    .tab-controls label{cursor:pointer; padding:.5rem .75rem; border:1px solid var(--border); border-radius:.5rem; background:#fff; font-weight:600; font-size:.95rem}
    /* hide tab radios but keep accessible */
    input[name="tab"]{position:absolute; left:-9999px}

    /* active tab styles */
    input#tab1:checked ~ .tab-controls label[for="tab1"],
    input#tab2:checked ~ .tab-controls label[for="tab2"],
    input#tab3:checked ~ .tab-controls label[for="tab3"]{
      background:var(--accent-weak); border-color:var(--accent); color:var(--accent)
    }

    .tab-panel{display:none;padding:1rem}
    input#tab1:checked ~ .tab-panels #panel1,
    input#tab2:checked ~ .tab-panels #panel2,
    input#tab3:checked ~ .tab-panels #panel3{display:block}

    /* =========================
       Cards & tables (generic)
       ========================= */
    .card{background:#fff;border:1px solid var(--border);border-radius:.75rem;padding:1rem}
    table{border-collapse:collapse;width:100%}
    th, td{border:1px solid var(--border);padding:.5rem;text-align:left;font-size:.95rem;vertical-align:top}
    th{background:#f3f4f6}
    .muted{color:var(--muted)}

    /* =========================
       Chips & length badges
       (reused in Sample tab)
       ========================= */
    .chips .chip{display:inline-block;margin:.15rem .3rem;padding:.15rem .45rem;border:1px solid var(--border);border-radius:.4rem;background:#f9fafb}
    .len-ok{background:#ecfdf5; border:1px solid #bbf7d0; padding:.1rem .3rem; border-radius:.3rem;}
    .len-below{background:#fef2f2; border:1px solid #fecaca; padding:.1rem .3rem; border-radius:.3rem;}
    .len-above{background:#fff7ed; border:1px solid #fed7aa; padding:.1rem .3rem; border-radius:.3rem;}

    /* =========================
       Panel 2 (no-JS host/geo)
       — radios + labels control spans inside the table you inject
       ========================= */
    /* visually hide functional radios (hostview + geoview) */
    .vh{
      position:absolute; width:1px; height:1px;
      margin:-1px; padding:0; border:0; overflow:hidden;
      clip:rect(0 0 0 0); clip-path:inset(50%); white-space:nowrap;
    }

    /* compact toolbar wrapper above table */
    #panel2 .sb-card{display:grid; grid-template-rows:auto 1fr; gap:.25rem}
    #panel2 .toolbars{display:flex; gap:.75rem; justify-content:flex-end; align-items:center; margin:0 0 .25rem 0; padding-right:.5rem}
    .controls{display:inline-flex; gap:.25rem; align-items:center}
    .controls .group-label{color:var(--muted); font-weight:600; font-size:.85rem; margin-right:.25rem}
    .controls label{
      display:inline-block; white-space:nowrap;
      padding:.15rem .5rem; font-size:.8rem; line-height:1;
      border:1px solid var(--border); border-radius:.4rem; background:#fff; cursor:pointer;
    }

    /* Hide all variants by default.
       Your injected table must use these span classes inside cells: 
       .host--species / .host--order / .host--class
       .geo--country / .geo--region / .geo--continent */
    #panel2 .sb-table-wrap .host{display:none}
    #panel2 .sb-table-wrap .geo{display:none}

    /* Show Host variant based on the checked radio (hostview) */
    #panel2 #hv-species:checked ~ .sb-card .sb-table-wrap .host--species{display:inline}
    #panel2 #hv-order:checked   ~ .sb-card .sb-table-wrap .host--order  {display:inline}
    #panel2 #hv-class:checked   ~ .sb-card .sb-table-wrap .host--class  {display:inline}

    /* Show Geography variant based on the checked radio (geoview) */
    #panel2 #gv-country:checked   ~ .sb-card .sb-table-wrap .geo--country  {display:inline}
    #panel2 #gv-region:checked    ~ .sb-card .sb-table-wrap .geo--region   {display:inline}
    #panel2 #gv-continent:checked ~ .sb-card .sb-table-wrap .geo--continent{display:inline}

    /* Highlight active labels (Host + Geography) */
    #panel2 #hv-species:checked ~ .sb-card .controls-host label[for="hv-species"],
    #panel2 #hv-order:checked   ~ .sb-card .controls-host label[for="hv-order"],
    #panel2 #hv-class:checked   ~ .sb-card .controls-host label[for="hv-class"],
    #panel2 #gv-country:checked   ~ .sb-card .controls-geo label[for="gv-country"],
    #panel2 #gv-region:checked    ~ .sb-card .controls-geo label[for="gv-region"],
    #panel2 #gv-continent:checked ~ .sb-card .controls-geo label[for="gv-continent"]{
      background:var(--accent-weak); color:var(--accent); border-color:var(--accent);
    }

    /* =========================
       Print: show all panels
       ========================= */
    @media print{
      .tab-controls{display:none}
      .tab-panel{display:block}
      .tabs{border:none}
    }
  </style>
</head>
<body>
  <header class="container">
    <h1><!--FILENAME--> Final Report</h1>
    <p class="muted"> HTML interactive report for Single-Sample mode</p>
  </header>

  <main class="container">
    <section class="tabs card">
      <!-- Tabs radios (no-JS switching) -->
      <input type="radio" id="tab1" name="tab" checked>
      <input type="radio" id="tab2" name="tab">
      <input type="radio" id="tab3" name="tab">

      <div class="tab-controls">
        <label for="tab1">Sample</label>
        <label for="tab2">Segment Background</label>
        <label for="tab3">Segment Mutations</label>
      </div>

      <div class="tab-panels">
        <!-- ============ TAB 1: SAMPLE ============ -->
        <article id="panel1" class="tab-panel">
          <h2>Sample</h2>
          <!--SAMPLE_CARD-->
          <!--SAMPLE_TABLE-->
        </article>

        <!-- ============ TAB 2: SEGMENT BACKGROUND ============ -->
        <article id="panel2" class="tab-panel">
          <h2>Segment Background</h2>

          <!-- Functional radios (hidden) that control Host & Geography views via CSS only -->
          <input class="vh" type="radio" id="hv-species" name="hostview" checked>
          <input class="vh" type="radio" id="hv-order"   name="hostview">
          <input class="vh" type="radio" id="hv-class"   name="hostview">

          <input class="vh" type="radio" id="gv-country"   name="geoview" checked>
          <input class="vh" type="radio" id="gv-region"    name="geoview">
          <input class="vh" type="radio" id="gv-continent" name="geoview">

          <div class="sb-card card">
            <!-- Host + Geography compact controls (labels only; radios above) -->
            <div class="toolbars">
              <div class="controls controls-host" aria-label="Host view">
                <span class="group-label">Host:</span>
                <label for="hv-species">Species</label>
                <label for="hv-order">Order</label>
                <label for="hv-class">Class</label>
              </div>
              <div class="controls controls-geo" aria-label="Geography view">
                <span class="group-label">Geography:</span>
                <label for="gv-country">Country</label>
                <label for="gv-region">Region</label>
                <label for="gv-continent">Continent</label>
              </div>
            </div>

            <!-- Inject here your full flat table (including the span variants inside cells) -->
            <div class="sb-table-wrap">
              <!--SEGMENT_BACKGROUND_TABLE-->
              <!--
                Expected structure inside the injected table (example cell content):
                  <td>
                    <span class="host host--species">{{HOST_SPECIES}}</span>
                    <span class="host host--order">{{HOST_ORDER}}</span>
                    <span class="host host--class">{{HOST_CLASS}}</span>
                  </td>
                  <td>
                    <span class="geo geo--country">{{GEO_COUNTRY}}</span>
                    <span class="geo geo--region">{{GEO_REGION}}</span>
                    <span class="geo geo--continent">{{GEO_CONTINENT}}</span>
                  </td>
                Length badges may use: <span class="len-ok|len-below|len-above">1234</span>
              -->
            </div>
          </div>
        </article>

        <!-- ============ TAB 3: SEGMENT MUTATIONS ============ -->
        <article id="panel3" class="tab-panel">
          <h2>Segment Mutations</h2>
          <!--SEGMENT_MUTATIONS_TABLE-->
        </article>
      </div>
    </section>
  </main>
</body>
</html>

"""

def to_html_report(sample_card, sample_table, seg_back_table, mut_table, html_skeleton,output_file_path,filename):
    """
    Generate the final HTML report by inserting the sample card and table into the HTML skeleton.
    """
    html_content = html_skeleton
    html_content = html_content.replace('<!--FILENAME-->', filename.replace("_final_report",""))
    html_content = html_content.replace('<!--SAMPLE_CARD-->', sample_card)
    html_content = html_content.replace('<!--SAMPLE_TABLE-->', sample_table)
    html_content = html_content.replace('<!--SEGMENT_BACKGROUND_TABLE-->', seg_back_table)
    html_content = html_content.replace('<!--SEGMENT_MUTATIONS_TABLE-->', mut_table)
    
    with open(os.path.join(output_file_path,f'{filename}.html'), 'w') as f:
        f.write(html_content)


def generate_final_report(df_segments, flags, seg_lens, mappings, muts_loci_meaning, html_skeleton, output_file_path,filename):
    """
    Generate the final HTML report using the provided data and flags.
    """
    sample_card = create_sample_card(flags)
    sample_table = create_sample_table(df_segments, flags, seg_lens, mappings)
    seg_back_table = create_segment_background_table(df_segments, flags, seg_lens, mappings)
    mut_table = create_segment_mutations_table(flags, muts_loci_meaning)
    
    to_html_report(sample_card, sample_table, seg_back_table, mut_table, html_skeleton,output_file_path,filename)