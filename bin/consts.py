TRNA_DICT = {
    'tRNA-Ala': 'trnA', 'tRNA-Cys':'trnC','tRNA-Asp':'trnD','tRNA-Glu':'trnE',
    'tRNA-Phe':'trnF','tRNA-Gly':'trnG','tRNA-His':'trnH','tRNA-Ile':'trnI',
    'tRNA-Lys':'trnK','tRNA-Leu':'trnL','tRNA-Met':'trnM','tRNA-Asn':'trnN',
    'tRNA-Pro':'trnP','tRNA-Gln':'trnQ','tRNA-Arg':'trnR','tRNA-Ser':'trnS',
    'tRNA-Thr':'trnT','tRNA-Val':'trnV','tRNA-Trp':'trnW','tRNA-Tyr':'trnY'}
REP_DICT={
          '1':'tRNA-Phe',
          '17':'tRNA-His',
          '14':'tRNA-Lys',
          '16':'tRNA-Arg',
          '13':'tRNA-Asp',
          '19':'tRNA-Glu',
          '12':'tRNA-Ser',
          '20':'tRNA-Thr',
          '9':'tRNA-Asn',
          '5':'tRNA-Gln',
          '8':'tRNA-Ala',
          '2':'tRNA-Val',
          '18':'tRNA-Leu',
          '11':'tRNA-Tyr',
          '7':'tRNA-Trp',
          '21':'tRNA-Pro',
          '15':'tRNA-Gly',
          '10':'tRNA-Cys',
          '6':'tRNA-Met',
          '4':'tRNA-Ile',
          '22':'12S ribosomal RNA',
          '23':'16S ribosomal RNA',
          '24':'ATP6',
          '25':'ATP8',
          '26':'COX1',
          '27':'COX2',
          '28':'COX3',
          '29':'CYTB',
          '30':'ND1',
          '31':'ND2',
          '32':'ND3',
          '33':'ND4',
          '34':'ND4L',
          '35':'ND5',
          '36':'ND6',
          '37':'ATP9',
          '38':'mutS',
          '39':'heg',
          '40':'secY',
          '41':'Reph',
          '42':'ORFX',
          '43':'RNAX'}
        
BASE_SET = {'trnA',
            'trnC',
            'trnD',
            'trnE',
            'trnF',
            'trnG',
            'trnH',
            'trnI',
            'trnK',
            'trnL',
            'trnL*',
            'trnM',
            'trnN',
            'trnP',
            'trnQ',
            'trnR',
            'trnS',
            'trnS*',
            'trnT',
            'trnV',
            'trnW',
            'trnY'}
# Import colormaps from plt
from mycolorpy import colorlist as mcp
TRNA_DICT = {
    'tRNA-Ala': 'trnA', 'tRNA-Cys':'trnC','tRNA-Asp':'trnD','tRNA-Glu':'trnE',
    'tRNA-Phe':'trnF','tRNA-Gly':'trnG','tRNA-His':'trnH','tRNA-Ile':'trnI',
    'tRNA-Lys':'trnK','tRNA-Leu':'trnL','tRNA-Met':'trnM','tRNA-Asn':'trnN',
    'tRNA-Pro':'trnP','tRNA-Gln':'trnQ','tRNA-Arg':'trnR','tRNA-Ser':'trnS',
    'tRNA-Thr':'trnT','tRNA-Val':'trnV','tRNA-Trp':'trnW','tRNA-Tyr':'trnY'}

REPLACEMENT_DICT = {'rrnS':['small subunit ribsomal RNA','small ribosomal RNA subunit RNA','12S small subunit ribosomal RNA','small ribosomal subunit RNA','12S large subunit ribosomal RNA''small ribosomal RNA subunit RNA','12S small subunit ribosomal RNA','small ribosomal RNA','small subunit rRNA','rRNA-12S','Small subunit ribosomal RNA','12S ribosomal RNA','12s ribosomal RNA','s-rRNA','rrnS','small subunit ribosomal RNA','12S rRNA','srRNA',
                                       '12S-rRNA','12S','ssu rRNA','rns','12srRNA','12s ribosomal RNA',
                                       'rnS','SSU','small ribosomal RNA subunit RNA', 'rrs', 'RNS', 'rns', 'RNS', 'rnsA', 'rrnaS', '12srRNA', '12SrRNA', 'rrn12S'],
                  'rrnL':['large subunit ribsomal RNA','large ribosomal RNA subunit RNA','16S large subunit ribosomal RNA','large ribosomal subunit RNA','16S large subunit ribosomal RNA''large ribosomal RNA subunit RNA','16S small subunit ribosomal RNA','large ribosomal RNA','large subunit rRNA','rRNA-16S','Large subunit ribosomal RNA','16S ribosomal RNA','16s ribosomal RNA','l-rRNA','rrnL','large subunit ribosomal RNA','16S rRNA',
                                       'lrRNA','16S-rRNA','16S','lsu rRNA','rnl','l6S ribosomal RNA',
                                       '16SrRNA','rrs','LSU','16S ribosomal RNA','rnL','large ribosonal RNA subunit RNA', 'rrl', 'RNL', 'rnl', 'rrnl', 'RNL', 'rnlA', 'rrnl', 'rnl_b','16srRNA', '16SrRNA', 'rrn16S'],
                  'cox1':['COX1','coI','cytochrome oxidase subunit-1','cytochrome c subunit I','cytochrome oxidase I','cytochrome C oxidase subunit 1','cytochrome oxidase c subunit 1','Cytochrome c oxidase subunit 1','cytochrome oxidase subunit I','cytochrome oxidase sununit 1','Cytochrome oxidase subunit I','cytochrome c oxidase subunit 1','COI','cytochrome c oxidase subunit I','CO1', 'cytochrome c oxidase subunit I COX1',
                          'cox1','COXI','coi','cytochrome oxidase subunit 1','cytochrome oxidase subunit I',
                          'Cox1','CO-I','coxI','COI protein'],
                  'cox2':['COX2','cytochrome oxidase subunit-2','cytochrome c subunit II','cytochrome oxidase II','cytochrome C oxidase subunit 2','cytochrome oxidase c subunit 2','Cytochrome c oxidase subunit 2','cytochrome oxidase subunit II','cytochrome oxidase sununit 2','Cytochrome oxidase subunit II','cytochrome c oxidase subunit 2','COII','cytochrome c oxidase subunit II','CO2', 'cytochrome c oxidase subunit II COX2',
                          'cox2','COXII','coii','cytochrome oxidase subunit 2','cytochrome oxidase subunit II',
                          'Cox2','CO-II','coxII','COII protein'],
                  'cox3':['COX3','CoxIII','cytochrome oxidase subunit-3','cytochrome c subunit III','cytochrome oxidase III','cytochrome C oxidase subunit 3','cytochrome oxidase c subunit 3','Cytochrome c oxidase subunit 3','cytochrome oxidase subunit III','cytochrome oxidase sununit 3','Cytochrome oxidase subunit III','cytochrome c oxidase subunit 3','COIII','cytochrome c oxidase subunit III','CO3', 'cytochrome c oxidase subunit IIII COX3',
                          'cox3','COXIII','coiii','cytochrome oxidase subunit 3','cytochrome c oxidase subunit 3',
                          'Cox3','CO-III','coxIII','COIII protein'],
                  'cob':['cytochrome b oxidase','cytochrome B','CYTB','apocytochrome b','cytochrome c oxidase subunit b','CytB','Cytochroome b','cob','cytb','cyt b','Cyt b',
                          'Cytb','Cyt B','cytochrome b','Cb','cyt-B','cytB','CYB','Cyt-b','Cytb protein'],
                  'nad1':['NADH1 dehydrogenase subunit 1','NADH dehydogenase subunit 1','ND1','NADH dehydrogenase subunit-1','cytochrome c oxidase subunits I','NADH dehydrogenase subunits 1','NADH dehydrogenase subunit1','NADH dehydrogenase subunit 1','NADH-ubiquinone oxidoreductase subunit 1','NADH dehyrogenase subunit 1','NADH dehydrogenase subunit I','NADH  dehydrogenase 1','NADH dehydrogenase 1','nad1-1','nad1-0','NADH dehydrogenase subunit 1','nad1','NADH1','nd1','nadh1','NAD1',
                         'NADH dehydrognase subunit I','nadh1','ndh1','Nad1','ND-1','ND1 protein','NADH hydrognase subunit 1','NDI', 'NADH-ubiquinone oxidoreductase chain 1'],
                  'nad2':['nad2_a','NADH2 dehydrogenase subunit 2','NADH dehydogenase subunit 2','ND2','NADH dehydrogenase subunit-2','cytochrome c oxidase subunits II','NADH dehydrogenase subunits 2','NADH dehydrogenase subunit2','NADH dehydrogenase subunit 2','NADH-ubiquinone oxidoreductase subunit 2','NADH dehyrogenase subunit 2','NADH dehydrogenase subunit II','NADH  dehydrogenase 2','NADH dehydrogenase 2','nad2-1','nad2-0','NADH dehydrogenase subunit 2','nad2','NADH2','nd2','nadh2','NAD2',
                         'NADH dehydrogenase subunit II','nadh2','ndh2','Nad2','ND-2','ND2 protein','NADH hydrogenase subunit 2','NDII', 'NADH-ubiquinone oxidoreductase chain 2'],
                  'nad3':['NADH3 dehydrogenase subunit 3','NADH dehydogenase subunit 3','ND3','NADH dehydrogenase subunit-3','cytochrome c oxidase subunits III','NADH dehydrogenase subunits 3','NADH dehydrogenase subunit3','NADH dehydrogenase subunit 3','NADH-ubiquinone oxidoreductase subunit 3','NADH dehyrogenase subunit 3','NADH dehydrogenase subunit III','NADH  dehydrogenase 3','NADH dehydrogenase 3','nad3-1','nad3-0','NADH dehydrogenase subunit 3','nad3','NADH3','nd3','nadh3','NAD3',
                         'NADH dehydrogenase subunit III','nadh3','ndh3','Nad3','ND-3','ND3 protein','NADH hydrogenase subunit 3','NDIII', 'NADH-ubiquinone oxidoreductase chain 3'],
                  'nad4':['NADH4 dehydrogenase subunit 4','NADH dehydogenase subunit 4','ND4','NADH dehydrogenase subunit-4','cytochrome c oxidase subunits IV','NADH dehydrogenase subunits 4','NADH dehydrogenase subunit4','NADH dehydrogenase subunit 4','NADH-ubiquinone oxidoreductase subunit 4','NADH dehyrogenase subunit 4','NADH dehydrogenase subunit IV','NADH  dehydrogenase 4','NADH dehydrogenase 4','nad4-1','nad4-0','NADH dehydrogenase subunit 4','nad4','NADH4','nd4','nadh4','NAD4',
                         'NADH dehydrogenase subunit IV','nadh4','ndh4','Nad4','ND-4','ND4 protein','NADH hydrogenase subunit 4','NDIV', 'NADH-ubiquinone oxidoreductase chain 4'],
                  'nad5':['NAHD dehydrogenase subunit 5','NADH5 dehydrogenase subunit 5','NADH dehydogenase subunit 5','ND5','NADH dehydrogenase subunit-5','cytochrome c oxidase subunits V','NADH dehydrogenase subunits 5','NADH dehydrogenase subunit5','NADH dehydrogenase subunit 5','NADH-ubiquinone oxidoreductase subunit 5','NADH dehyrogenase subunit 5','NADH dehydrogenase subunit V','NADH  dehydrogenase 5','NADH dehydrogenase 5','nad5-1','nad5-0','NADH dehydrogenase subunit 5','nad5','NADH5','nd5','nadh5','NAD5',
                         'NADH dehydrogenase subunit V','nadh5','ndh5','Nad5','ND-5','ND5 protein','NADH hydrogenase subunit 5','NDV', 'NADH-ubiquinone oxidoreductase chain 5'],
                  'nad6':['NADH6 dehydrogenase subunit 6','NADH dehydogenase subunit 6','ND6','NADH dehydrogenase subunit-6','cytochrome c oxidase subunits VI','NADH dehydrogenase subunits 6','NADH dehydrogenase subunit6','NADH dehydrogenase subunit 6','NADH-ubiquinone oxidoreductase subunit 6','NADH dehyrogenase subunit 6','NADH dehydrogenase subunit VI','NADH  dehydrogenase 6','NADH dehydrogenase 6','nad6-1','nad6-0','NADH dehydrogenase subunit 6','nad6','NADH6','nd6','nadh6','NAD6',
                         'NADH dehydrogenase subunit VI','nadh6','ndh6','Nad6','ND-6','ND6 protein','NADH hydrogenase subunit 6','NDVI', 'NADH-ubiquinone oxidoreductase chain 6'],
                  'atp8':['ATP synthetase F0 subunit 8','ATP synthasee subunit 8','ATP8','ATPase subunit-8','ATPase subunits 8','ATPase subunit-8','ATPase subunits 8','ATP synthase 8','ATP synthase FO subunit 8','ATP synthase F0 subunit 8','atp8','ATPase8','ATPase 8', 'ATP synthase subunit 8',
                          'ATP synthase subunit 8','ATP 8','AT8','atpase8','ATP synthetase subunit 8','adenosine triphosphate subunit 8','Atp8','ATPase subunit 8','ATPase sunuint 8'],
                  'atp6':['ATP synthetase F0 subunit 6','ATP synthasee subunit 6','ATP6','ATPase subunit-6','ATPase subunits 6','ATPase subunit-6','ATPase subunits 6','ATP synthase 6','ATP synthase FO subunit 6','ATP synthase F0 subunit 6','atp6','ATPase6','ATPase 6', 'ATP synthase subunit A','ATP 6','atpase6','ATP synthetase subunit 6','adenosine triphosphate subunit 6','Atp6','ATPase subunit 6','ATPase subuint 6', 'ATP synthase subunit alpha'
                          'ATP synthase subunit 6','ATP 6','AT6','atpase6','ATP synthetase subunit 6','adenosine triphosphate subunit 6','Atp6','ATPase subunit 6','ATPase subuint 6'],
                  'nad4L':['NADH4L dehydrogenase subunit 4L','NADH dehydogenase subunit 4L','ND4L','NADH dehydrogenase subunit-4L','cytochrome c oxidase subunits IV L','NADH dehydrogenase subunits 4L','NADH dehydrogenase subunit4L','NADH dehydrogenase subunit 4 L','NADH-ubiquinone oxidoreductase subunit 4L','NADH dehyrogenase subunit 4L','NADH dehydrogenase subunit IV L','NADH  dehydrogenase 4L','NADH dehydrogenase 4L','NADH dehydrogenase subunit 4L','nad4l','nad4L','NADH4L','nd4L','NADH dehydrognase subunit 4L', 'NADH-ubiquinone oxidoreductase chain 4L',
                          'nadh4L','NAD4L','NADH dehydrognase subunit IV L','nadh4l','Nad4L','ND-4L','ND4L protein','NADH hydrogenase subunit 4L','ND4l','NDIVL','NADH dehydrogenase subunit 4l'],
                  'atp9':['ATP9','ATPase subunits 9','ATPase subunit-9','ATPase subunits 9','ATP synthase 9','ATP synthase F0 subunit 9','atp9','ATPase9','ATPase 9','ATP synthase subunit 9','ATP 9',
                          'AT9','atpase9','ATP synthetase subunit 9','adenosine triphosphate subunit 9','Atp9',
                          'ATPase subunit 9','ATPase subuint 9'],
                  'mutS':['msh1','DNA mismatch repair protein MutS','MutS-like protein','putative mismatch repair protein','DNA mismatch repair protein mutS',
                          'DNA mismatch repair protein','MutS-like protein','mismatch repair protein'],
                  'heg':['homing endonuclease','truncated homing endonuclease','homing endonuclease'],
                  'nqo1' : ['NADH-quinone oxidoreductase'],
                  'secY':['SecY-independent transporter protein'],
                  'Reph':['replication helicase subunit'],
                  'ORFX':[],
                  'nad7':['ND7', 'NADH dehydrogenase subunit 7'],
                  'nad9':['ND9', 'nad9_2', 'NADH dehydrogenase subunit 9'],
                  'RNAX':[],
                  'nad8':['ND8'],
                  'nad10':['ND10', 'NADH dehydrogenase subunit 10'],
                  'nad11':['ND11', 'NADH dehydrogenase subunit 11'],
                  'cox11':['COX11'],
                  'atp1':['ATP1'],
                  'atp3':['ATP3'],
                  'rrp1' : ['RRP1', 'ribosomal protein L1'],
                  'rrp2' : ['RRP2', 'ribosomal protein L2'],
                  'rrp3' : ['RRP3', 'ribosomal protein L3'],
                  'rrp4' : ['RRP4', 'ribosomal protein L4'],
                  'rrp5' : ['RRP5', 'ribosomal protein L5'],
                  'rrp6' : ['RRP6', 'ribosomal protein L6'],
                  'rrp7' : ['RRP7', 'ribosomal protein L7'],
                  'rrp8' : ['RRP8', 'ribosomal protein L8'],
                  'rrp9' : ['RRP9', 'ribosomal protein L9'],
                  'rrp10' : ['RRP10', 'ribosomal protein L10'],
                  'rrp11' : ['RRP11', 'ribosomal protein L11'],
                  'rrp12' : ['RRP12', 'ribosomal protein L12'],
                  'rrp13' : ['RRP13', 'ribosomal protein L13'],
                  'rrp14' : ['RRP14', 'ribosomal protein L14'],
                  'rrp15' : ['RRP15', 'ribosomal protein L15'],
                  'rrp16' : ['RRP16', 'ribosomal protein L16'],
                  'rrp17' : ['RRP17', 'ribosomal protein L17'],
                  'rrp18' : ['RRP18', 'ribosomal protein L18'],
                  'rrs1': ['RRS1', 'ribosomal protein S1'],
                  'rrs2': ['RRS2', 'ribosomal protein S2'],
                  'rrs3': ['RRS3', 'ribosomal protein S3'],
                  'rrs4': ['RRS4', 'ribosomal protein S4'],
                  'rrs5': ['RRS5', 'ribosomal protein S5'],
                  'rrs6': ['RRS6', 'ribosomal protein S6'],
                  'rrs7': ['RRS7', 'ribosomal protein S7'],
                  'rrs8': ['RRS8', 'ribosomal protein S8'],
                  'rrs9': ['RRS9', 'ribosomal protein S9'],
                  'rrs10': ['RRS10', 'ribosomal protein S10'],
                  'rrs11': ['RRS11', 'ribosomal protein S11'],
                  'rrs12':['RRS12', 'ribosomal protein S12'],
                  'rrs13':['RRS13', 'ribosomal protein S13'],
                  'rrs14':['RRS14', 'ribosomal protein S14'],
                  'rrs15':['RRS15', 'ribosomal protein S15'],
                  'rrs16':['RRS16', 'ribosomal protein S16'],
                  'rrs17':['RRS17', 'ribosomal protein S17'],
                  'rrs18':['RRS18', 'ribosomal protein S18'],}
MTDNA_GENES = ['nad1', 'nad2', 'nad3', 'nad4', 'nad6', 'nad4L', 'nad5', 'nad8', 'nad9', 'nad11','cox1', 'cox2', 'cox3', 'cob', 'rnl', 'rns', 'atp8','atp9', 'atp6', 'rrnL', 'rrnS', 'RNR1', 'RNR2'] # list of 'classic' mtDNA genes

COLORS = {'Extreme (32C)': '#e74c3c', 'Normal (25C)': '#f1c40f', 'Low (15C)': '#2ecc71'}

MTDNA_GENES_TO_COMPLEX = {
'I' : ['nad1', 'nad2', 'nad3', 'nad4', 'nad5', 'nad6', 'nad4L', 'nad8', 'nad9', 'nad10', 'nad11'],
'III' : ['cox1', 'cox2', 'cox3', 'cox4', 'cox5', 'cox6', 'cox7', 'cox8', 'cox9'],
'IV' : ['cob'],
'V' : ['atp8', 'atp9', 'atp6']}
MTDNA_GENES_TO_COMPLEX_FLAT = {}
for key in MTDNA_GENES_TO_COMPLEX:
    for i in MTDNA_GENES_TO_COMPLEX[key]:
        MTDNA_GENES_TO_COMPLEX_FLAT[i] = key

STRAND_NAMES = {'+':True, '-':False, '1':True, '2':False, 'F':True, 'R':False, 'forward':True, 'reverse':False, 'f':True, 'r':False, '1' : True, '-1' : False, 1 : True, -1 : False, 0 : False, '0' : False}

DM_TRANSCRIPTS = [
  [19405, 6405, True],
  [9530, 12134, True],
  [3716, 93, False],
  [10891, 6349, False],
  [15077, 11894, False]
]

DM = 'Drosophila melanogaster'

RANKS = ['domain','kingdom','phylum','class','order','family','genus']

atp6_atp8_hairpin_hs = 'AAUGAACGAAAAUCUGUUCGCUUCAUU'
atp6_atp8_hairpin_hs_structure = '((((((((((......)))).))))))'

SYNS = {
            'N' : ['A','C','G','T'],
            'R' : ['G','A'],
            'Y' : ['T','C'],
            'M' : ['A','C'],
            'K' : ['G','T'],
            'S' : ['G','C'],
            'W' : ['A','T'],
            'H' : ['A','T','C'],
            'B' : ['T','G','C'],
            'D' : ['T','G','A'],
            'V' : ['A','G','C'],
            'A' : ['A'],
            'C' : ['C'],
            'G' : ['G'],
            'T' : ['T']
            }
CODON_TABLES={2:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S",                           # noqa: E241
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "AGA":"STOP","AGG":"STOP","TAA":"STOP","TAG":"STOP"
    },
    5:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
        4:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
        9:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
        13:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "G", "AGG": "G",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
        14:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Y",               # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAG":"STOP"},
        21:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "N", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "S",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
       24:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                           # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "K",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAA":"STOP","TAG":"STOP"},
        33:{
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "Y",               # noqa: E241
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "S", "AGG": "K",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "TAG":"STOP"
        }}
PHYLUM_TO_CODON_TABLE = {
        'Chordata': 2,      # Vertebrate Mitochondrial
        'Arthropoda': 5,    # Invertebrate Mitochondrial
        'Echinodermata': 9, # Echinoderm Mitochondrial
        'Mollusca': 5,      # Invertebrate Mitochondrial
        'Nematoda': 5,      # Invertebrate Mitochondrial
        'Annelida': 5,      # Invertebrate Mitochondrial
        'Platyhelminthes': 9, # Echinoderm Mitochondrial
        'Cnidaria': 4,      # Mold Mitochondrial
        'Porifera': 4,      # Mold Mitochondrial
        'Ascomycota': 4,    # Mold Mitochondrial
        'Basidiomycota': 4, # Mold Mitochondrial
    }

INVALID_SPECIES ="""Oxyurichthys formosanus
Amblema plicata
Paraglypturus tonganus
Urosaurus nigricaudus
Corydoras hastatus
Vestiaria coccinea
Chlorodrepanis stejnegeri
Ctenopharyngodon idellus
Mesoclemmys hogei
Uragus sibiricus
Anodonta lucida
Anodonta arcaeformis
Gloydius saxatilis
Tylorrhynchus heterochaetus
Glossogobius giuris
Ceratopsyche cerva
Ceratopsyche trifora
Acanthocobitis zonalternans
Erynnis montanus
Sinotaia purificata
Alopecoenas salamonis
Sinapis arvensis
Corydoras trilineatus
Garrulax perspicillatus
Metrioptera bonneti
Cricetulus migratorius
(Cyprinus carpio
Melophus lathami
Tanakia tanago
Lamprologus kungweensis
Corydoras cruziensis
Candida auris
Stichaeus grigorjewi
Candida ipomoeae
Gallirallus striatus
Diploprion bifasciatum
Ceratopsyche gautamittra
Tribolodon brandtii
Tribolodon sachalinensis
Lamprologus meleagris
Oreoglanis macropterus
Paracentropogon rubripinnis
Pomatomus saltatrix
Kazachstania sinensis
Amazilia versicolor
Parophasma galinieri
Loliolus uyii
Microphysogobio yaluensis
Etheostoma chuckwachatte
Bathyclupea megaceps
Grus leucogeranus
Garrulax formosus
Saprochaete suaveolens
Saprochaete ingens
Kazachstania unispora
Tachycineta meyeni
Doryichthys boaja
Asymblepharus himalayanus
Fasciolopsis buski
Ceratopsyche columnata
Taeniura meyeni
Hoplobatrachus chinensis
Hyporhamphus limbatus
Ganoderma flexipes
Saprochaete fungicola
Sahyadria chalakkudiensis
Porzana pusilla
Candida vartiovaarae
Candida norvegica
Paraneetroplus synspilus
Stachybotrys chlorohalonata
Thranita danae
Solenaia oleivora
Nitzschia sp.
Esanthelphusa keyini
Odontobutis potamophila
Synagrops philippinensis
Ficus variegata
Garrulax chinensis
Hylocharis cyanus
Garrulax albogularis
Barbodes lateristriga
Nemachilichthys rueppelli
Catostomus discobolus
Porzana paykullii
Sepia aculeata
Pennahia macrocephalus
Elaphe taeniurus
Etheostoma camurum
Etheostoma tippecanoe
Amaurornis akool
Echinosciurus variegatoides
Syntheosciurus granatensis
Cyprinus multitaeniata
Padda oryzivora
Hemigrammus bleheri
Diploderma micangshanensis
Cryodraco antarcticus
Hadrosciurus ignitus
Acosmetura nigrogeniculata
Lechriodus melanopyga
Candida sake
Spinibarbus denticulatus
Corydoras duplicareus
Corydoras arcuatus
Corydoras panda
Ophiocara porocephala
Amazilia brevirostris
Corydoras sterbai
Humphaplotropis culaishanensis
Phymatostetha huangshanensis
Cricetulus kamensis
Anadara sativa
Cobitis matsubarai
Candida psychrophila
Megalaima virens
Tariqlabeo bicornis
Paradactylodon mustersi
Candida intermedia
Gloydius brevicaudus
Proedromys liangshanensis
Mesobuthus martensii
Eothenomys chinensis
Leucoraja erinacea
Triplophysa stoliczkai
Ambystoma unisexual
Tribolodon nakamurai
Ctenoptilum vasava
Leptobotia mantschurica
Candida alai
Aceros waldeni
Crassostrea gigas
Chrysochir aureus
Sylvia crassirostris
Crassostrea angulata
Microtus kikuchii
Sardinops melanostictus
Nakaseomyces castellii
Nakaseomyces bacillisporus
Crocodylus johnsoni
Microtus fortis
Schistosoma spindale
Ustilago maydis
Kazachstania servazzii
Corydoras rabauti
Ceratopsyche fukienensis
Diplectrona hexapetala
Kisaura zhejiangensis
Candida railenensis
Sinularia flexibilis
Microtus montebelli
Locustella pleskei
Aphanius iberus
Corydoras aeneus
Corydoras paleatus
Neolamprologus similis
Candida subhashii
Hyphessobrycon anisitsi
Corydoras pygmaeus
Garrulax courtoisi
Hebius metusium
Barilius ardens
Niwaella nigrolinea
Microdous chalmersi
Oreonectes anophthalmus
Troglonectes daqikongensis
Troglonectes furcocaudalis
Troglonectes jiarongensis
Troglonectes shuilongensis
Margaritifera homsensis
Corydoras concolor
Corydoras julii
Microtus maximowiczii
Myodes rutilus
Peltigera neocanina
Thamnaconus multilineatus
Equetus lanceolatus
Ariomma melanum
Rosa hybrid
Daedalea dickinsii
Fulgoraria rupestris
Canarium labiatum
Eothenomys olitor
Fissocantharis imparicornis
Phyllothelys sinensis
Asymblepharus sikimmensis
Anchoa tricolor
Orthopristis ruber
Ferrissia fragilis"""
INVALID_SPECIES = INVALID_SPECIES.split('\n')

RCRS_ACCS = 'NC_012920'

RCRS_SEQ = """GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT
ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCA
AACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAA
ACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC
TTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATA
CCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAA
GCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTC
CTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGT
TCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTC
AAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAA
ACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGC
GGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCC
TCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGAC
TACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGA
TACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAA
CACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGG
AGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATA
CCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAG
ACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAG
AAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTAAG
AGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTC
AAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGT
CGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACA
AAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTA
GCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAA
AGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATG
AAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAA
TTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCT
ACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATA
GGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAG
TTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTC
CAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTA
ACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCA
CTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATC
ACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGC
CTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAAC
AAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAA
AAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGC
ATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCT
AACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCC
ACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCG
GGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTA
CCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGA
CCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAA
CTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACA
GCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCA
GGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTAC
GTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACNTTCAAATTCCTCC
CTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGA
TATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTTGTTAAGATGGC
AGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTT
AACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAATCGCAATGGCA
TTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAAC
GTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAA
GAGCCCCTAAAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCT
CTCACCATCGCTCTTCTACTATGAACCCCCCTCCCCATACCCAACCCCCTGGTCAACCTC
AACCTAGGCCTCCTATTTATTCTAGCCACCTCTAGCCTAGCCGTTTACTCAATCCTCTGA
TCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGCGCACTGCGAGCAGTAGCCCAA
ACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGC
TCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCA
TGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTC
GACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGC
CCCTTCGCCCTATTCTTCATAGCCGAATACACAAACATTATTATAATAAACACCCTCACC
ACTACAATCTTCCTAGGAACAACATATGACGCACTCTCCCCTGAACTCTACACAACATAT
TTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCC
CGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTA
GCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAACC
TAAGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCC
CCTTATTTCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCC
ACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCG
AAAATGTTGGTTATACCCTTCCCGTACTAATTAATCCCCTGGCCCAACCCGTCATCTACT
CTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAG
TAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAAAATAAACCCTC
GTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTC
TAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCA
ATCAATACTCATCATTAATAATCATAATAGCTATAGCAATAAAACTAGGAATAGCCCCCT
TTCACTTCTGAGTCCCAGAGGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTC
TCACATGACAAAAACTAGCCCCCATCTCAATCATATACCAAATCTCTCCCTCACTAAACG
TAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAGCAGGCAGTTGAGGTGGATTAA
ACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAA
TAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATCC
TAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTAT
CTCGCACCTGAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCC
TAGGAGGCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCA
CAAAAAACAATAGCCTCATCATCCCCACCATCATAGCCACCATCACCCTCCTTAACCTCT
ACTTCTACCTACGCCTAATCTACTCCACCTCAATCACACTACTCCCCATATCTAACAACG
TAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCCCACACTCATCG
CCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTTA
GGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGT
AACAGCTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAA
TTAAGCTAAGCCCTTACTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCT
AAGCACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAG
AAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTC
GGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCACTCA
GCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACA
AAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTC
TAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACA
TCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCA
TCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGG
CGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGC
TCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAG
CAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAG
CAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAAC
CCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAG
TCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACC
GCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAAC
ACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAA
TAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTA
TGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATAT
TTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCG
CTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGA
AATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCC
TGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACG
TTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCT
TCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCC
ATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCC
TATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCC
TATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTT
GAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGG
AGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAAT
CTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCC
TCCATGACTTTTTCAAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAAT
TATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACGC
TACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTT
CCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAA
TACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCAT
CCTAGTCCTCATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGA
TCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTACTGAACCTACGAGTACACCGA
CTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCAGGCGA
CCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTAT
AATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAAC
AGATGCAATTCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATA
CTACGGTCAATGCTCTGAAATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGA
ATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTA
CCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAG
AGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTATGGCCCACCAT
AATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAA
CTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGA
ACCAAAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCC
GCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATC
AACAACCGACTAATCACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATA
ACCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTT
ATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTA
TCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCACAGTGATTATAGGCTTTCGC
TCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATC
CCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTA
CGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACC
CTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTA
CTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTA
AGCCTCTACCTGCACGACAACACATAATGACCCACCAATCACATGCCTATCATATAGTAA
AACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAG
CCATGTGATTTCACTTCCACTCCATAACGCTCCTCATACTAGGCCTACTAACCAACACAC
TAACCATATACCAATGATGGCGCGATGTAACACGAGAAAGCACATACCAAGGCCACCACA
CACCACCTGTCCAAAAAGGCCTTCGATACGGGATAATCCTATTTATTACCTCAGAAGTTT
TTTTCTTCGCAGGATTTTTCTGAGCCTTTTACCACTCCAGCCTAGCCCCTACCCCCCAAT
TAGGAGGGCACTGGCCCCCAACAGGCATCACCCCGCTAAATCCCCTAGAAGTCCCACTCC
TAAACACATCCGTATTACTCGCATCAGGAGTATCAATCACCTGAGCTCACCATAGTCTAA
TAGAAAACAACCGAAACCAAATAATTCAAGCACTGCTTATTACAATTTTACTGGGTCTCT
ATTTTACCCTCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCA
TCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGGACTTCACGTCATTATTGGCT
CAACTTTCCTCACTATCTGCTTCATCCGCCAACTAATATTTCACTTTACATCCAAACATC
ACTTTGGCTTCGAAGCCGCCGCCTGATACTGGCATTTTGTAGATGTGGTTTGACTATTTC
TGTATGTCTCCATCTATTGATGAGGGTCTTACTCTTTTAGTATAAATAGTACCGTTAACT
TCCAATTAACTAGTTTTGACAACATTCAAAAAAGAGTAATAAACTTCGCCTTAATTTTAA
TAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACTACCACAACTCA
ACGGCTACATAGAAAAATCCACCCCTTACGAGTGCGGCTTCGACCCTATATCCCCCGCCC
GCGTCCCTTTCTCCATAAAATTCTTCTTAGTAGCTATTACCTTCTTATTATTTGATCTAG
AAATTGCCCTCCTTTTACCCCTACCATGAGCCCTACAAACAACTAACCTGCCACTAATAG
TTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTGACTAC
AAAAAGGATTAGACTGAACCGAATTGGTATATAGTTTAAACAAAACGAATGATTTCGACT
CATTAAATTATGATAATCATATTTACCAAATGCCCCTCATTTACATAAATATTATACTAG
CATTTACCATCTCACTTCTAGGAATACTAGTATATCGCTCACACCTCATATCCTCCCTAC
TATGCCTAGAAGGAATAATACTATCGCTGTTCATTATAGCTACTCTCATAACCCTCAACA
CCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAG
CAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGTAC
ATAACCTAAACCTACTCCAATGCTAAAACTAATCGTCCCAACAATTATATTACTACCACT
GACATGACTTTCCAAAAAACACATAATTTGAATCAACACAACCACCCACAGCCTAATTAT
TAGCATCATCCCTCTACTATTTTTTAACCAAATCAACAACAACCTATTTAGCTGTTCCCC
AACCTTTTCCTCCGACCCCCTAACAACCCCCCTCCTAATACTAACTACCTGACTCCTACC
CCTCACAATCATGGCAAGCCAACGCCACTTATCCAGTGAACCACTATCACGAAAAAAACT
CTACCTCTCTATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGA
ACTAATCATATTTTATATCTTCTTCGAAACCACACTTATCCCCACCTTGGCTATCATCAC
CCGATGAGGCAACCAGCCAGAACGCCTGAACGCAGGCACATACTTCCTATTCTACACCCT
AGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACT
AAACATTCTACTACTCACTCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTT
AATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTT
ATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGT
ACTCTTAAAACTAGGCGGCTATGGTATAATACGCCTCACACTCATTCTCAACCCCCTGAC
AAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTC
CATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACAT
AGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCAT
TCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTA
CGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACT
AATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAA
CCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCT
ACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAAC
ACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAA
CACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCAT
TACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAA
CAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCC
CATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAG
GCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACC
ACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAA
AAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTC
TTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGA
GCCACAACCCAAACAACCCAGCTCTCCCTAAGCTTCAAACTAGACTACTTCTCCATAATA
TTCATCCCTGTAGCATTGTTCGTTACATGGTCCATCATAGAATTCTCACTGTGATATATA
AACTCAGACCCAAACATTAATCAGTTCTTCAAATATCTACTCATCTTCCTAATTACCATA
CTAATCTTAGTTACCGCTAACAACCTATTCCAACTGTTCATCGGCTGAGAGGGCGTAGGA
ATTATATCCTTCTTGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCC
ATTCAAGCAATCCTATACAACCGTATCGGCGATATCGGTTTCATCCTCGCCTTAGCATGA
TTTATCCTACACTCCAACTCATGAGACCCACAACAAATAGCCCTTCTAAACGCTAATCCA
AGCCTCACCCCACTACTAGGCCTCCTCCTAGCAGCAGCAGGCAAATCAGCCCAATTAGGT
CTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCAC
TCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAA
AATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCA
GCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGT
CAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCAC
ATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCAC
AACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTC
ACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGT
TTCTACTCCAAAGACCACATCATCGAAACCGCAAACATATCATACACAAACGCCTGAGCC
CTATCTATTACTCTCATCGCTACCTCCCTGACAAGCGCCTATAGCACTCGAATAATTCTT
CTCACCCTAACAGGTCAACCTCGCTTCCCCACCCTTACTAACATTAACGAAAATAACCCC
ACCCTACTAAACCCCATTAAACGCCTGGCAGCCGGAAGCCTATTCGCAGGATTTCTCATT
ACTAACAACATTTCCCCCGCATCCCCCTTCCAAACAACAATCCCCCTCTACCTAAAACTC
ACAGCCCTCGCTGTCACTTTCCTAGGACTTCTAACAGCCCTAGACCTCAACTACCTAACC
AACAAACTTAAAATAAAATCCCCACTATGCACATTTTATTTCTCCAACATACTCGGATTC
TACCCTAGCATCACACACCGCACAATCCCCTATCTAGGCCTTCTTACGAGCCAAAACCTG
CCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAG
CACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTC
CTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAACCTATTCCCCCG
AGCAATCTCAATTACAATATATACACCAACAAACAATGTTCAACCAGTAACTACTACTAA
TCAACGCCCATAATCATACAAAGCCCCCGCACCAATAGGATCCTCCCGAATCAACCCTGA
CCCCTCTCCTTCATAAATTATTCAGCTTCCTACACTATTAAAGTTTACCACAACCACCAC
CCCATCATACTCTTTCACCCACAGCACCAATCCTACCTCCATCGCTAACCCCACTAAAAC
ACTCACCAAGACCTCAACCCCTGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCGC
TGTAGTATATCCAAAGACAACCATCATTCCCCCTAAATAAATTAAAAAAACTATTAAACC
CATATAACCTCCCCCAAAATTCAGAATAATAACACACCCGACCACACCGCTAACAATCAA
TACTAAACCCCCATAAATAGGAGAAGGCTTAGAAGAAAACCCCACAAACCCCATTACTAA
ACCCACACTCAACAGAAACAAAGCATACATCATTATTCTCGCACGGACTACAACCACGAC
CAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACG
CAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATC
CAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAAT
CACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATC
AATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAA
TGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGG
ATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGC
AACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGT
AATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTG
AGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTT
GCCCTTCATTATTGCAGCCCTAGCAACACTCCACCTCCTATTCTTGCACGAAACGGGATC
AAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACAC
AATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATT
CTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCC
TCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCC
TAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCC
CATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTA
TTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTA
CCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCT
AATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCTGTCCTTGTAGTAT
AAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGA
GAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCT
CTGTTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACA
ACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTACCATAAAT
ACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTA
CAAGCAAGTACAGCAATCAACCCTCAACTATCACACATCAACTGCAACTCCAAAGCCACC
CCTCACCCACTAGGATACCAACAAACCTACCCACCCTTAACAGTACATAGTACATAAAGC
CATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCC
TCAGATAGGGGTCCCTTGACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCT
ACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACAT
CTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGAC
ATCACGATG""".replace("\n", "")

ND4_start_position_minus1 = 10_759
CYTB_start_position_minus1 = 14746
COX1_start_position_minus1 = 5903

GENE_START_DICT = {
        'ND4' : ND4_start_position_minus1 + 1,
        'CYTB' : CYTB_start_position_minus1 + 1,
        'COX1' : COX1_start_position_minus1 + 1,
        'RNR2' : 1671,
        'RNR1' : 684,
        'TS2_TL2_ND5' : 12207
        
}