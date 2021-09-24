#!/share/apps/python/bin/python
import sys, os

workingdir = os.getcwd()
#datadir = workingdir.split('cafe')[0] + '/cafe/src/data/'
datadir = workingdir + '/data/'

# fasta
hm_fasta = datadir + 'fasta/human/hg19/hg19.fa'
ms_fasta = datadir + 'fasta/mouse/mm9/mm9.fa'

# flag table
flag_table = datadir + 'info/flag_table'

# matrix
hm_matrix = datadir + 'matrix/human/'
ms_matrix = datadir + 'matrix/mouse/'

# cage
hm_fpseqAll = datadir + 'anno/cage/human/hg19/hg19.cage_peak_coord_robust.bed'
ms_fpseqAll = datadir + 'anno/cage/mouse/mm9/mm9.cage_peak_coord_robust.bed'

# polya
hm_tpseqAll = datadir + 'anno/polya/human/hg19/hg19_all.15.sumCM.bed'
ms_tpseqAll = datadir + 'anno/polya/mouse/mm9/mm9_all.15.sumCM.bed'
