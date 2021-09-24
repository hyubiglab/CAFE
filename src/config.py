#!/share/apps/python/bin/python
import sys, os

# genome
hm_chrs = map(lambda x: 'chr' + x, map(str, range(1, 23)) + ['X', 'Y'])
hm_intronmin = 61; hm_intronmax = 265006

ms_chrs = map(lambda x: 'chr' + x, map(str, range(1, 20)) + ['X', 'Y'])
ms_intronmin = 52; ms_intronmax = 240764

# sampling
sampling_number = 1000000 #sampling number

# rpds
rpds_type = 'hmm' #rpds type ('dis' or 'hmm')
rpds_order = 3 #rpds order (1~5)

# exj
exj_type = 'signal' #exj type ('read' or 'signal')
exj_cutoff = '0.217' #exj cutoff

# tss
tss_type = 'all' #tss type ('all', 'major', or 'near')

hm_fpseqInterval = int(3000)
hm_fpseqThreshold = float(0.05)

ms_fpseqInterval = int(2000)
ms_fpseqThreshold = float(0.05)

# cps
cps_type = 'all' #cps type ('all', major', or 'near')

hm_tpseqInterval = int(5000)
hm_tpseqThreshold = float(0.05)

ms_tpseqInterval = int(3000)
ms_tpseqThreshold = float(0.05)
