#!/share/apps/python/bin/python
import sys, os

import config as conf
import module as module

assembly = sys.argv[1]
gtfFile = sys.argv[2]
bamFile = sys.argv[3]
outputdir = sys.argv[4]

if assembly == 'hg19': chrs = conf.hm_chrs ##
elif assembly == 'mm9': chrs = conf.ms_chrs ##

for chr in chrs:
	gtf = module.getGtf(gtfFile)
	trxs = gtf[chr]
	outputFile = open(outputdir + '/ratio_' + chr + '.txt', 'w')
	for i in xrange(len(trxs)):
		trx = trxs[i]; exons = trx.exons()
		for exon in exons:
			strandRatio, senseReads, antisenseReads = module.checkRatioExon(bamFile, exon)
			outputFile.write(trx.geneid() + '\t' + trx.trxid() + '\t' + exon.coord() + '\t' + exon.sense() + '\t' + str(len(senseReads)) + '\t' + str(len(antisenseReads)) + '\t' + str(strandRatio) + '\n')
	outputFile.close()
