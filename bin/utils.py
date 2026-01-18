from cmath import inf
import pandas as pd
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez
import os, re
import snoop
import consts
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from importlib import reload
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu 
from statannotations.Annotator import Annotator
import mygene

reload(consts)
PATH = os.getcwd()
Entrez.email = 'your_email@example.com'
Entrez.API = '[REMOVED_API_KEY]' 

# Initialize once
mg = mygene.MyGeneInfo()

def enst_to_symbol(enst_id):
    """Convert ENSEMBL transcript ID to gene symbol using mygene."""
    try:
        result = mg.query(enst_id, 
                         scopes='ensembl.transcript',
                         fields='symbol',
                         species='human')
        if result['hits']:
            return result['hits'][0].get('symbol', pd.NA)
        return enst_id
    except:
        return enst_id

def enst_to_symbol_batch(transcript_ids):
    """Vectorized conversion using mygene batch query."""
    results = mg.querymany(transcript_ids.tolist(),
                          scopes='ensembl.transcript',
                          fields='symbol',
                          species='human',
                          returnall=True)
    
    # Create mapping dict
    id_to_symbol = {}
    for res in results['out']:
        enst = res['query']
        symbol = res.get('symbol', pd.NA)
        if symbol is pd.NA:
            id_to_symbol[enst] = enst
        else:
            id_to_symbol[enst] = symbol
    
    return transcript_ids.map(id_to_symbol)

def gene_symbol_to_number(to_rep,number2symbol = False):
    """ recives a str and dict, replaces str with any matching values in dict and returns dict reversed to replace string annotations
    with numerical annotations
    
    Parameters
    ----------
    to_rep : str
        string to replace
    number2symbol : bool
        if true, returns dict reversed
    
    Returns
    -------
    str
        string with replaced values
    dict
        dict with replaced values
    """
    counter=0
    bools=False
    if number2symbol:
        temp_rep = {v:k for k,v in consts.REP_DICT.items()}
    else:
        temp_rep = consts.REP_DICT
    if '*' in to_rep: #checks if this gene is a duplicate marked by *
        counter=to_rep.count('*')
    to_rep=to_rep.replace('*','')
    if '-' in to_rep[0]: #checks if the gene is in the complementary strand
        bools=True
        to_rep=to_rep.replace('-','',1)
    for k,v in temp_rep.items(): #iterates over the replacement dictionary and looks for matching gene strings
        if to_rep==v:
            to_rep=k
            break
    to_rep=to_rep+'*'*counter
    if bools: to_rep='-'+to_rep #If gene is in complementry, if the answer is yes 
    return to_rep

def list_of_genes_to_number(genes, number2symbol = False):
    return [gene_symbol_to_number(i, number2symbol = number2symbol) for i in genes]

def row_iter_to_df(df, func, *args, **kwargs):
    """
    Iterates over rows of a dataframe and applies a function to each row.
    Returns a dataframe with the results.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to iterate over.
    func : function
        Function to apply to each row.
    *args : list
        Arguments to pass to func.
    **kwargs : dict 
        Keyword arguments to pass to func.
    
    Returns
    -------
    pandas.DataFrame
        Dataframe with results of func applied to each row.
    
    """
    total_df = pd.DataFrame()
    for _, row in df.iterrows():
        temp = pd.DataFrame(func(row, *args, **kwargs))
        if type(temp) != pd.DataFrame: continue
        if total_df.size == 0: total_df = temp
        else: total_df = total_df.append(temp)
    return total_df

def dir_checker(dirname):
  """
  check if a directory exists, if it does not, create it.

  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  if not os.path.isdir(os.path.join(PATH, dirname)):
    os.mkdir(os.path.join(PATH, dirname))

def reverse_complement(seq):
    """
    """
    return str(Seq.Seq(seq).reverse_complement())

def move_legend(ax, new_loc, **kws):
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, loc=new_loc, title=title, **kws)

def record_check(ID, type = 'gb'):
  """
  Check whether org genbank exists, if not download it.
  """
  suffix = '.gbk' if type == 'gb' else '.fasta'
  loc = 'genbank_DB' if type == 'gb' else 'fasta_DB'
  filename = os.path.join(PATH, loc, ID + suffix)
  if not os.path.isfile(filename):
    print('Downloading ' + ID + '...')
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = type, retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()
    print('Done!')

def genomic_ranges(sample : pd.DataFrame, no_trna = False) -> list:
    """
    Based on a sample df made by the annotation.py script, return a list in the following fromat:
    [start, end, gene_name]
    Parameters
    ----------
    sample : pd.DataFrame
      A dataframe made using the annotation.py code. MUST HAVE COLUMNS:
      Gene, Position
    not_trna : bool
      True if junctions ignore tRNAs 
    Returns
    -------
    granges : list
      List of all the granges in the following format:
      [start, end, gene_name]
    """
    no_na = sample.Gene.unique().tolist()
    try: no_na.remove(np.nan) #remove nan if it exists
    except ValueError: pass 
    mtdna_len = len(sample) # get the length of the mtDNA
    if no_trna:
      genelist = [gene for gene in no_na if 'trn' not in gene.lower()]     
    else:
      genelist = no_na
    granges = []
    for gene in genelist: # for each gene in the gene list
        cur_gene = sample.loc[sample.Gene == gene, :]
        if len(cur_gene) < 50: continue
        cur_pos = cur_gene.Position.to_list()
        zero_ind = cur_pos[0]
        last_ind = cur_pos[-1]
        size = abs(last_ind - zero_ind)
        if  size > (mtdna_len* 0.9): # This indicates that the gene is too long to be considered because it wraps around the mtDNA
          print(size)
          last_ind =  max([i for i in cur_pos if i < (mtdna_len* 0.9)])
          print(last_ind)
        granges.append([zero_ind, last_ind, gene])
    return granges
  
def wrap(start, end, pos, flank, size):
    """
    If the current position + the flanking region or the current position - the flanking region exceeds the 
    boundaries of the mtDNA, wrap around (because the mtDNA is circular)

    Parameters
    ----------
    start : int
      Gene start position
    end : int
      Gene ending position
    pos : int
      Current position (in the loop)
    flank : int
      The size of the flanking region we want
    size : int
      The mtDNA genome size
    
    Returns 
    -------
    bool
      A boolean that is True if the current position wraps around the mtDNA and is flanking.
    """
    if start - flank < 0:
        return pos >= start - flank + size
    elif start + flank > size:
        return pos < start + flank - size
    elif end - flank < 0:
        return pos >= end - flank + size
    elif end + flank > size:
        return pos < end + flank - size
    else: return False

def downstream_wrap(end, pos, flank, size):
  if end + flank > size:
    return pos < end + flank - size   
  else: return False

def flanking(pos, genelist, flank, gsize, direction = 'both'):
    """
    Receive a genelist, return True if the current position is within flanking range (including wraps)
    
    Parameters
    ----------
    pos : int
      The current position to check
    genelist : list
      The list to compare the position with, MUST be in the following format: [start, end, gene_name]
    flank : int
      The flanking region size we are interested in.
    gsize : int
      The size of the mtDNA genome.
    
    Returns
    -------
    bool
      True if the pos is within flanking range.
    """
    if direction == 'both':
      return any([(grange[0] - flank <= pos < grange[0] + flank) or (grange[1] - flank <= pos < grange[1] + flank) or (wrap(grange[0], grange[1], pos, flank, gsize)) for grange in genelist])
    elif direction == 'downstream':
      return any([(grange[1] < pos < grange[1] + flank) or (downstream_wrap(grange[1], pos, flank, gsize)) for grange in genelist])
    elif direction == 'upstream':pass #TODO(Noam): add upstream option.

def Which_tRNA(name, verbosity = 1):
    """ This code recieves a gene name, tries to match it with common tRNA motifs, and then with any matching values in the replacement dict, if nothing matches, returns None
    which is later removed"""
    name = name.replace('*', '')
    if name in consts.REPLACEMENT_DICT.keys(): return name 
    if 'trn' in name: return name
    c=0
    name=str(name)
    if 'transfer' in name:
        name = name.replace('transfer ','t')
    codon = re.search(r'\(\w{3}\)',name) # Grab codon if it exists
    if codon: codon = codon.group(0)
    else: codon=''
    name = re.sub(pattern=r'\(\w{3}\)',repl = '',string = name)
    if name in consts.TRNA_DICT.keys():  # if the tRNA is in tRNA-X format, replace with my format and return
        name = consts.TRNA_DICT[name]
        return name+codon

    #Following lines try to match with each possible tRNA.
    if re.match(r't[r,R]\D{1,2}(H|His)',name):
        name='tRNA-His'
        
    if re.match(r't[r,R]\D{1,2}(K|Lys)',name):
        name='tRNA-Lys'
        
    if re.match(r't[r,R]\D{1,2}(R|Arg)',name):
        name='tRNA-Arg'
        
    if re.match(r't[r,R]\D{1,2}(D|Asp)',name):
        name='tRNA-Asp'
        
    if re.match(r't[r,R]\D{1,2}(E|Glu)',name):
        name='tRNA-Glu'
        
    if re.match(r't[r,R]\D{1,2}S',name):
        name='tRNA-Ser'
        
    if re.match(r't[r,R]\D{1,2}T',name):
        name='tRNA-Thr'
        
    if re.match(r't[r,R]\D{1,2}(N|Asn)',name):
        name='tRNA-Asn'
        
    if re.match(r't[r,R]\D{1,2}(Q|Gln)',name):
        name='tRNA-Gln'
        
    if re.match(r't[r,R]\D{1,2}V',name):
        name='tRNA-Val'
        
    if re.match(r't[r,R]\D{1,2}L',name):
        name='tRNA-Leu'
        
    if re.match(r't[r,R]\D{1,2}I',name):
        name='tRNA-Ile'
        
    if re.match(r't[r,R]\D{1,2}M',name):
        name='tRNA-Met'
        
    if re.match(r't[r,R]\D{1,2}(F|Phe)',name):
        name='tRNA-Phe'
        
    if re.match(r't[r,R]\D{1,2}(Y|Tyr)',name):
        name='tRNA-Tyr'
        
    if re.match(r't[r,R]\D{1,2}(W|Trp)',name):
        name='tRNA-Trp'
        
    if re.match(r't[r,R]\D{1,2}P',name):
        name='tRNA-Pro'
        
    if re.match(r't[r,R]\D{1,2}G',name):
        name='tRNA-Gly'
        
    if re.match(r't[r,R]\D{1,2}C',name):
        name='tRNA-Cys'
    
    if re.match(r't[r,R]\D{1,2}A',name):
        name='tRNA-Ala'
        
    for k,v in consts.REPLACEMENT_DICT.items(): # if not tRNA, this loop matches according to replacement dict
        if name in v:
            name=k
            return name
        if c==len(consts.REPLACEMENT_DICT)-1:
            name=None
            return name
        c=+1
    try:
        if name in consts.TRNA_DICT.keys(): # if matched with any of the if's above (for tRNAs), replace with my annotation accordingly.
            name = consts.TRNA_DICT[name]
            return name+codon
    except KeyError:
        pass
    if verbosity > 0:
        print(name)
    if 'RNA' in name or 'rna' in name:
        return None
    elif 'orf' in name or 'ORF' in name:
        return 'ORFX'
    else:
        if verbosity > 0:
            print(f"Could not replace {name} - returning it")
        return name

    
def mtdna_region(pos, window, total, left, inclusive = True) -> str:
    """
    Designed for circular DNA, return a list of positions (INDEX 1 BASED) to the left or right side of pos

    Parameters
    ----------
    pos : int
        The current position to return a window around
    window : int
        The size of the window
    total : int
        The total size of the DNA
    left : bool
        True if the window should be to the left of the
    inclusive : bool
        True if the current position should be included in the window
        (Default value = True)

    Returns
    -------
    positions : list
        A list of positions (INDEX 1 BASED) to the left or right side of pos
    """

    if type(left) != bool: raise TypeError(f'Parameter left must be either True (left side window) or False (right side window! Given {left} instead')
    if total < window or total < pos:
        raise ValueError(f'The total size of the DNA must be smaller than both the window and the position!\nParameters given:\nTotal = {total}\nPosition = {pos}\nWindow = {window}')
    
    positions = []
    if left:
        if pos - window < 1:
            positions += list(range(total - abs(window - pos) + + (1 if inclusive else 0), total + 1))
            positions += list(range(1, pos + (1 if inclusive else 0)))
        else:
            positions += list(range(pos - window + (1 if inclusive else 0), pos + (1 if inclusive else 0)))
    else: #right
        if pos + window > total:
            positions += list(range(1, (pos + window) - total + (0 if inclusive else 1)))
            positions += list(range(pos + (0 if inclusive else 1), total + 1))
        else:
            positions += list(range(pos + (0 if inclusive else 1), pos + window + 1))
    return positions

def checkInt(str):
    if type(str) != str: return str
    if str[0] in ('-', '+'):
        return str[1:].isdigit()
    return str.isdigit()

def remove_codon(gene):
    """
    Remove the codon from the gene name
    """
    return re.sub(pattern='\(\w{3}\)',repl = '',string = gene)

def conv_gene_loc_format(loc):
    """
    Convert the gene locations to the format: [start, end, strand] from a string of start:end:strand
    """
    loc = loc.split(':')
    if len(loc) <= 2:
        loc = ['0', '0', '1']   
    loc = [0 if not checkInt(j) else int(re.sub(string = j, pattern = '\D+', repl = '')) for j in loc] 
    return loc

def remove_trna(gorder, glocs):
    """
    Remove the tRNA from the gene name
    """
    new_gorder = []
    new_glocs = []
    for i, gene in enumerate(gorder):
        if 'trn' not in gene:
            new_glocs.append(glocs[i])
            new_gorder.append(gene)
    return new_gorder, new_glocs

def detect_overlapping_feature(start, end, strand, features, feature_locations, stranded_mode = False):
    """
    Detect if a feature overlaps with a given locus. If it does, return the feature name. If it doesn't, return a string of the two closest features

    Parameters
    ----------
    start : int
        The start position of the feature
    end : int
        The end position of the feature
    features : list
        A list of features to check against
    feature_locations : dict
        A list of the location of features, organized in the same order as the features list in the format: [start:end:strand]
    stranded_mode : bool
        True if the strand of the feature must be the same as the strand of the locus. False if the strand of the feature doesn't matter.

    Returns
    -------
    feature : str
        The name of the feature that overlaps with the current feature
    """
    features = [remove_codon(gene) for gene in features]
    if strand in consts.STRAND_NAMES:
        strand = consts.STRAND_NAMES[strand]
    feature = None
    start_dists = []
    end_dists = []
    for i, f in enumerate(features):
        loc = feature_locations[i]
        if type(loc) != list:
            loc = conv_gene_loc_format(loc)
            feature_locations[i] = loc
        f_strand = loc[2]
        if f_strand in consts.STRAND_NAMES:
            f_strand = consts.STRAND_NAMES[f_strand]
        if start >= loc[0] and end <= loc[1] and (strand == f_strand if stranded_mode else True):
            # If the feature starts before the current feature ends and ends after the current feature starts, and the strand is the same as the current feature (in stranded mode), then the current feature overlaps with the current feature
            feature = f
            return feature
        if i == 0:
            end_dist = 50000
        if i == len(features)-1:
            start_dist = 50000
        start_dist = start - loc[0] if start > loc[0] else np.inf
        end_dist = loc[1] - end if end < loc[1] else np.inf
        start_dists.append(start_dist)
        end_dists.append(end_dist)
    left_gene = features[np.argmin(start_dists)]
    right_gene = features[np.argmin(end_dists)]
    return f'{left_gene}_{right_gene}'

def wrapper_for_detect_overlapping_feature(row, stranded_mode = False, no_trna = True):
    start = row['start']
    end = row['stop']
    strand = row['strand']
    if no_trna:
        features = row['gorder_notrna']
        feature_locations = row['gloc_notrna']
    else:
        features = row['Gene_order']
        feature_locations = row['Gene_locations']
    return detect_overlapping_feature(start, end, strand, features, feature_locations, stranded_mode)

def generate_motif_list(start, end):
    """
    Generate a list of motifs from a given start and end position

    Parameters
    ----------
    start : int
        The start position of the motif
    end : int
        The end position of the motif
    
    Returns
    -------
    motif_list : list
        A list of motifs from start to end
    
    """
    start = int(start)
    end = int(end)
    region = mtdna_region(start, window = end - start, total = 100000, left = False, inclusive = True)
    return region

def generate_motif_list_wrapper(row):
    """
    Wrapper for generate_motif_list
    """
    start = row['start']
    end = row['stop']
    return generate_motif_list(start, end)

def filter_for_coding_regions(record):
    """
    
    """
    coding_seq = ''
    seq = record.seq
    for feature in record.features:
        try:
            if feature.type == 'CDS':
                cur_seq = feature.extract(seq)
                if feature.strand == -1:
                    cur_seq = cur_seq.reverse_complement()
                coding_seq += cur_seq
        except:
            print(f'Error extracting coding sequence {record.id}')
            return ''
    return coding_seq

def record_load(ID, type = 'genbank'):
    """
    Load a record from the genbank database
    """
    if type == 'genbank':
        record = SeqIO.read(os.path.join(PATH, 'genbank_DB', f'{ID}.gbk'), 'genbank')
    else:
        record = SeqIO.read(os.path.join(PATH, 'fasta_DB', f'{ID}.fasta'), 'fasta')
    return record

def get_aa_sequence(ID, gene):
    """
    Recieve a RefSeq ID, extract the AA sequence of a given protein coding gene
    """
    record_check(ID)
    record = SeqIO.read(os.path.join(PATH, 'genbank_DB', f'{ID}.gbk'), 'genbank')
    try:
        for feat in record.features:
            if feat.type == 'CDS':
                if feat.qualifiers['gene'][0] == gene:
                    return feat.qualifiers['translation'][0]
    except KeyError:
        print(f'{gene} not found in {ID}')
        return ''

def plot_y_per_x_compare_cat_per_grp(
    df, x, y, cat, grps,
    grp_col = 'phylum',
    count_col = 'organism',
    ylab = 'RSCU',
    xlab = 'Amino acid',
    savefig = False):
    """
    Plot y per x for each category in cat, for each group in grps and compare the groups using a statistical test (Mann-Whitney U test) generates a plot inline

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing the data
    x : str
        The name of the column containing the x values
    y : str
        The name of the column containing the y values
    cat : str
        The name of the column containing the categories
    grps : list
        A list of the groups to compare
    grp_col : str, optional
        The name of the column containing the groups, by default 'phylum'
    count_col : str, optional
        The name of the column containing the counts, by default 'organism'
    ylab : str, optional
        The label for the y axis, by default 'Amino acid'
    xlab : str, optional
        The label for the x axis, by default 'RSCU'
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes object
    """
    for ph in grps:
        try:
            _, ax = plt.subplots(figsize=(6,4))
            data = df[df[grp_col] == ph]
            # If data is empty, skip to the next phylum
            if len(data) == 0:
                continue
            pairs = [[(i, True), (i, False)] for i in data[x].unique()]
            AAs = sorted(data[x].unique())
            sns.barplot(ax = ax, x = x, y = y, data = data, hue = cat, palette = 'tab10', order = AAs)
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.set_title(f'{ph} N = {len(data[count_col].unique())}')
            annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=AAs, hue = cat, verbose=False)
            annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', comparisons_correction = 'benjamini-hochberg', correction_format = 'replace', fontsize = 8)
            annotator.apply_and_annotate() 
            # Move legend
            move_legend(ax, "upper left") 
            plt.tight_layout()

            # Despine the top and right borders
            sns.despine(ax = ax, offset = 10, trim = True)
            
            cur_pvals = multipletests([mannwhitneyu(data[(data.In_mtDNA == True) & (data['AA'] == i)][y], data[(data.In_mtDNA == False) & (data['AA'] == i)][y]).pvalue for i in AAs if\
                len(data[(data.In_mtDNA == True) & (data['AA'] == i)]) > 0 and len(data[(data.In_mtDNA == False) & (data['AA'] == i)]) > 0], method = 'fdr_bh')[1]  # For each amino acid within current phylum, perform a mann-whitney u test and get the pvalue and then perform a FDR correction.
            
            if savefig:
                plt.savefig(os.path.join(PATH, 'figures', f'supp_fig_s2_rscu_AA_phylum_{ph}_{savefig}.svg'), pad_inches = .7, dpi = 300)
        except ValueError:
            print(f'No data for {ph}')
            print('Skipping...')
        plt.show()
    
def parse_vienna_output(file_path):
    """
    Parse the output of the ViennaRNA package
    """
    atp6_atp8_Large_dict = {}
    with open(file_path, 'r') as f:
        atp6_atp8_Large = f.readlines()
        for line in atp6_atp8_Large:
            line = line.rstrip('\n')
            if line.startswith('>'):
                header = line.replace('>', '').split('_')[0] + '_' + line.replace('>', '').split('_')[1]
                atp6_atp8_Large_dict[header] = []
            else:
                if not re.search(r'[a-zA-Z]', line):
                    structure = line.split(' ')[0]
                    delta_g = line.split(' ')[1].replace('(', '').replace(')', '')
                    atp6_atp8_Large_dict[header].append(structure)
                    atp6_atp8_Large_dict[header].append(delta_g)
                else:
                    atp6_atp8_Large_dict[header].append(line)
        
        atp6_atp8_Large_df = pd.DataFrame.from_dict(atp6_atp8_Large_dict, orient='index', columns = ['overlap_seq', 'structure', 'delta_g']).reset_index().rename(columns={'index': 'species'})
        # Drop delta_g values of ''
        atp6_atp8_Large_df = atp6_atp8_Large_df[atp6_atp8_Large_df['delta_g'] != '']
        atp6_atp8_Large_df['delta_g'] = atp6_atp8_Large_df['delta_g'].astype(float)
        return atp6_atp8_Large_df

def add_polya(seq, expected_length=None):
    """
    Add poly(A) tail to complete a TAA stop codon for mitochondrial genes.

    Some mitochondrial genes have incomplete stop codons that are completed
    post-transcriptionally. This function adds 1 or 2 'A's as needed.

    If a TAA stop codon cannot be generated (likely due to SNP/indel in the
    last codon), the function adds the number of A's required to match the
    expected length (determined from the reference sequence).

    Parameters:
    -----------
    seq : str
        DNA sequence
    expected_length : int, optional
        Expected length of the sequence after adding poly(A) tail.
        This is determined from the reference sequence (NC_012920.1).
        If None, will add 2 A's when TAA cannot be generated.

    Returns:
    --------
    str : Sequence with poly(A) tail added (if needed)

    Notes:
    ------
    - If sequence already ends with TAA: returns as-is (no A's added)
    - If sequence ends with TA: adds 1 A to complete TAA
    - If sequence ends with T and adding 2 A's creates TAA: adds 2 A's
    - If TAA cannot be generated and expected_length is provided:
      adds the number of A's to match expected_length
    - If TAA cannot be generated and expected_length is None:
      adds 2 A's for consistency (assumes SNP/indel in last codon)
    """
    seq = seq.upper()  # Ensure uppercase for consistency

    # Check if already ends with TAA stop codon
    if seq.endswith('TAA'):
        return seq

    # Check if ends with TA (needs 1 A to complete TAA)
    if seq.endswith('TA'):
        seq_modified = seq + 'A'
        return seq_modified

    # Check if ends with T (potentially needs 2 A's to complete TAA)
    if seq.endswith('T') and not seq.endswith('TA'):
        seq_modified = seq + 'AA'
        if seq_modified.endswith('TAA'):
            return seq_modified
        else:
            # Could not create TAA, use expected length if provided
            if expected_length is not None:
                num_as_needed = expected_length - len(seq)
                seq_modified = seq + ('A' * num_as_needed)
                print(
                    f"WARNING: Could not generate TAA stop codon (sequence ends with T). "
                    f"Added {num_as_needed} A's to match expected length from reference. "
                    f"Sequence ends with: {seq_modified[-3:]}"
                )
                return seq_modified
            else:
                # Fallback: add 2 A's for consistency
                seq_modified = seq + 'AA'
                print(
                    f"WARNING: Could not generate TAA stop codon (sequence ends with T). "
                    f"Added 2 A's for consistency (assuming SNP/indel in last codon). "
                    f"Sequence ends with: {seq_modified[-3:]}"
                )
                return seq_modified

    # For sequences not ending with T, TA, or TAA, try adding A's
    seq_with_one_a = seq + 'A'
    if seq_with_one_a.endswith('TAA'):
        return seq_with_one_a

    seq_with_two_as = seq + 'AA'
    if seq_with_two_as.endswith('TAA'):
        return seq_with_two_as

    # If adding 2 A's doesn't create TAA, use expected length if provided
    if expected_length is not None:
        num_as_needed = expected_length - len(seq)
        seq_modified = seq + ('A' * num_as_needed)
        print(
            f"WARNING: Could not generate TAA stop codon. "
            f"Sequence ends with: {seq[-3:] if len(seq) >= 3 else seq}. "
            f"Added {num_as_needed} A's to match expected length from reference. "
            f"Modified sequence ends with: {seq_modified[-3:] if len(seq_modified) >= 3 else seq_modified}"
        )
        return seq_modified
    else:
        # Fallback: add 2 A's for consistency
        seq_modified = seq + 'AA'
        print(
            f"WARNING: Could not generate TAA stop codon. "
            f"Sequence ends with: {seq[-3:] if len(seq) >= 3 else seq}. "
            f"Added 2 A's for consistency (assuming SNP/indel in last codon). "
            f"Modified sequence ends with: {seq_modified[-3:] if len(seq_modified) >= 3 else seq_modified}"
        )
        return seq_modified


def fasta_to_df(fasta_path, polya = False, reference_id = 'NC_012920.1'):
    """
    Convert a fasta file to a dataframe - add poly A when necessary

    Parameters:
    -----------
    fasta_path : str
        Path to the FASTA file
    polya : bool, optional
        If True, add poly(A) tails to sequences
    reference_id : str, optional
        ID of the reference sequence (default: 'NC_012920.1')
        Used to determine expected sequence length when adding poly(A) tails

    Returns:
    --------
    pd.DataFrame : DataFrame with columns 'ID' and 'sequence'
    """
    with open(fasta_path, 'r') as f:
        fasta = f.readlines()
    fasta_dict = {}
    for line in fasta:
        line = line.rstrip('\n')
        if line.startswith('>'):
            header = line.split(' ')[0].replace('>', '')
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line
    fasta_df = pd.DataFrame.from_dict(fasta_dict, orient='index', columns = ['sequence']).reset_index().rename(columns={'index': 'ID'})

    if polya:
        # First, determine expected length from reference sequence
        expected_length = None
        if reference_id in fasta_df['ID'].values:
            reference_seq = fasta_df[fasta_df['ID'] == reference_id]['sequence'].iloc[0]
            # Apply add_polya to reference to determine expected length
            reference_with_polya = add_polya(reference_seq, expected_length=None)
            expected_length = len(reference_with_polya)
            print(f"\nReference sequence ({reference_id}) expected length: {expected_length}")
        else:
            print(f"\nWARNING: Reference sequence '{reference_id}' not found in FASTA file.")
            print("Will use fallback logic (add 2 A's when TAA cannot be generated).")

        # Apply add_polya to all sequences using the expected length
        fasta_df['sequence'] = fasta_df['sequence'].apply(lambda seq: add_polya(seq, expected_length=expected_length))
    return fasta_df

def count_non_standard_starts(fasta_path):
    """
    Counts sequences in a FASTA file that start with a character 
    other than A, G, C, or T (case-insensitive).
    """
    count = 0
    with open(fasta_path, 'r') as f:
        check_next_line = False
        
        for line in f:
            line = line.strip()
            if not line: continue # Skip empty lines
            
            if line.startswith(">"):
                # We found a header, so the next valid line is the start of the sequence
                check_next_line = True
            elif check_next_line:
                # This is the first line of sequence data
                first_char = line[0].upper()
                if first_char not in "AGCT":
                    count += 1
                
                # We have checked the start of this sequence; 
                # stop checking until we hit the next header '>'
                check_next_line = False

    return count


def gene_titv_to_dataframe(gene_titv_results):
    """
    Convert gene_titv_results dictionary to a single dataframe for supplementary tables.

    Parameters
    ----------
    gene_titv_results : dict
        Dictionary containing gene-level Ti/Tv analysis results.

    Returns
    -------
    pd.DataFrame
        Single dataframe with per-gene rows and a 'Total' summary row.
    """
    gene_results = gene_titv_results['gene_results']
    global_stats = gene_titv_results['global_gene_stats']

    rows = []
    for gene_name, gene_data in gene_results.items():
        row = {
            'Gene': gene_name,
            'Length (bp)': gene_data['length'],
            'N Codons': gene_data['n_codons'],
            'N Sequences': gene_data['n_sequences'],
            'Pos1 Ti': gene_data['first_pos']['ti_count'],
            'Pos1 Tv': gene_data['first_pos']['tv_count'],
            'Pos1 Ti/Tv': round(gene_data['first_pos']['ti_tv_ratio'], 2),
            'Pos1 π': gene_data['first_pos']['nucleotide_diversity'],
            'Pos2 Ti': gene_data['second_pos']['ti_count'],
            'Pos2 Tv': gene_data['second_pos']['tv_count'],
            'Pos2 Ti/Tv': round(gene_data['second_pos']['ti_tv_ratio'], 2),
            'Pos2 π': gene_data['second_pos']['nucleotide_diversity'],
            'Pos3 Ti': gene_data['third_pos']['ti_count'],
            'Pos3 Tv': gene_data['third_pos']['tv_count'],
            'Pos3 Ti/Tv': round(gene_data['third_pos']['ti_tv_ratio'], 2),
            'Pos3 π': gene_data['third_pos']['nucleotide_diversity'],
            'Pos1+2 Ti/Tv': round(gene_data['combined_12_titv'], 2),
            'Pos1+2 π': gene_data['combined_12_diversity'],
            'Ti/Tv Diff (Pos3 - Pos1+2)': round(gene_data['titv_diff'], 2),
            'Shows Protein Pattern': gene_data['shows_protein_pattern']
        }
        rows.append(row)

    # Add total row
    rows.append({
        'Gene': 'Total',
        'Length (bp)': sum(gene_results[g]['length'] for g in gene_results),
        'N Codons': sum(gene_results[g]['n_codons'] for g in gene_results),
        'N Sequences': '-',
        'Pos1 Ti': global_stats['first_pos_ti'],
        'Pos1 Tv': global_stats['first_pos_tv'],
        'Pos1 Ti/Tv': round(global_stats['first_pos_titv'], 2),
        'Pos1 π': global_stats['roa_diversity_first'],
        'Pos2 Ti': global_stats['second_pos_ti'],
        'Pos2 Tv': global_stats['second_pos_tv'],
        'Pos2 Ti/Tv': round(global_stats['second_pos_titv'], 2),
        'Pos2 π': global_stats['roa_diversity_second'],
        'Pos3 Ti': global_stats['third_pos_ti'],
        'Pos3 Tv': global_stats['third_pos_tv'],
        'Pos3 Ti/Tv': round(global_stats['third_pos_titv'], 2),
        'Pos3 π': global_stats['roa_diversity_third'],
        'Pos1+2 Ti/Tv': round(global_stats['combined_12_titv'], 2),
        'Pos1+2 π': global_stats['roa_diversity_12'],
        'Ti/Tv Diff (Pos3 - Pos1+2)': round(global_stats['third_pos_titv'] - global_stats['combined_12_titv'], 2),
        'Shows Protein Pattern': global_stats['shows_protein_pattern']
    })

    return pd.DataFrame(rows)


def combined_results_to_dataframe(combined_results):
    """
    Convert combined_results dictionary (MDP analysis) to a single dataframe.

    Parameters
    ----------
    combined_results : dict
        Dictionary containing MDP Ti/Tv analysis results.

    Returns
    -------
    pd.DataFrame
        Single dataframe with all MDPs from both rRNAs, including comparison statistics.
    """
    rows = []

    for rrna_key in ['rnr1_results', 'rnr2_results']:
        if rrna_key not in combined_results:
            continue

        rrna_data = combined_results[rrna_key]
        rrna_name = 'RNR1' if 'rnr1' in rrna_key else 'RNR2'
        mdp_results = rrna_data['mdp_results']
        comparison = rrna_data['comparison']

        for mdp_name, mdp_data in mdp_results.items():
            comp_data = comparison.get(mdp_name, {})
            row = {
                'MDP': mdp_name,
                'rRNA': rrna_name,
                'Start': mdp_data['start'],
                'End': mdp_data['end'],
                'Length (bp)': mdp_data['length'],
                'N Sequences': mdp_data['n_sequences'],
                'Pos1 Ti': mdp_data['first_pos']['ti_count'],
                'Pos1 Tv': mdp_data['first_pos']['tv_count'],
                'Pos1 Ti/Tv': round(mdp_data['first_pos']['ti_tv_ratio'], 2),
                'Pos1 π': mdp_data['first_pos']['nucleotide_diversity'],
                'Pos2 Ti': mdp_data['second_pos']['ti_count'],
                'Pos2 Tv': mdp_data['second_pos']['tv_count'],
                'Pos2 Ti/Tv': round(mdp_data['second_pos']['ti_tv_ratio'], 2),
                'Pos2 π': mdp_data['second_pos']['nucleotide_diversity'],
                'Pos3 Ti': mdp_data['third_pos']['ti_count'],
                'Pos3 Tv': mdp_data['third_pos']['tv_count'],
                'Pos3 Ti/Tv': round(mdp_data['third_pos']['ti_tv_ratio'], 2),
                'Pos3 π': mdp_data['third_pos']['nucleotide_diversity'],
                'Pos1+2 Ti/Tv': round(mdp_data['combined_12_titv'], 2),
                'Pos1+2 π': mdp_data['combined_12_diversity'],
                'Ti/Tv Diff (Pos3 - Pos1+2)': round(comp_data.get('mdp_titv_diff', mdp_data['third_pos']['ti_tv_ratio'] - mdp_data['combined_12_titv']), 2),
                'Random Mean Ti/Tv Diff': round(comp_data.get('random_mean_titv_diff', 0), 2) if comp_data else None,
                'Random Protein Pattern %': round(comp_data.get('random_protein_pattern_proportion', 0) * 100, 1) if comp_data else None,
                'p-value': comp_data.get('p_value') if comp_data else None,
                'Significant': comp_data.get('significant') if comp_data else None,
                'Shows Protein Pattern': comp_data.get('mdp_shows_protein_pattern', mdp_data['third_pos']['ti_tv_ratio'] > mdp_data['combined_12_titv'])
            }
            rows.append(row)

    return pd.DataFrame(rows)