#!/share/apps/python/bin/python
import sys, os

import module as module

assembly = sys.argv[1]
gtfFile = sys.argv[2]
chr = sys.argv[3]
outputdir = sys.argv[4]

def filterIncluTrxs(trxs):
	ntrxs = []
	trxNum = len(trxs)
	for i in xrange(trxNum):
		trx = trxs[i]; exons = trx.exons(); introns = trx.introns()
		if trx.trxid().find('TSS') > 0 and trx.trxid().find('CPS'):
			ntrxs.append(trx)
		else:
			ans = False
			trxToCheck = trxs[max(0, i-50):i]
			trxToCheck += trxs[i+1:min(trxNum, i+50)]
			for trx2 in trxToCheck:
				exons2 = trx2.exons(); introns2 = trx2.introns()
				# the parent trx should have one more exon to avoid messup
				if len(exons2) >= len(exons):
					prt = False
					if trx.trxid().find('TSS') > 0:
						if trx.sense() == '+':
							if trx.start() == trx2.start(): prt = True
						else:
							if trx.end() == trx2.end(): prt = True
					if trx.trxid().find('CPS') > 0:
						if trx.sense() == '+':
							if trx.end() == trx2.end(): prt = True
						else:
							if trx.start() == trx2.start(): prt = True
					if trx.trxid().find('TSS') < 0 and trx.trxid().find('CPS') < 0:
						prt = True
					if prt:
						inc = 0
						for exon in exons:
							for exon2 in exons2:
								if exon2.inclusive(exon): inc += 1
						inc2 = 0
						for intron in introns:
							for intron2 in introns2:
								if intron.start() == intron2.start() and intron.end() == intron2.end(): inc2 += 1
						if inc == len(exons) or inc2 == len(introns): ans = True; break
			if not ans: ntrxs.append(trx)

	return ntrxs

gtfFile = module.getGtf(gtfFile)
trxs = gtfFile[chr]; trxs.sort(module.cmp0)

ftrxs = []
trxNum = len(trxs)
for i in xrange(trxNum):
	trx = trxs[i]; exons = trx.exons()
	if len(exons) == 1:
#		if trx.trxid().find('TSS') > 0 or trx.trxid().find('CPS') > 0: ##
		if trx.trxid().find('TSS'): ##
			ftrxs.append(trx)
	else:
		ftrxs.append(trx)

ftrxs.sort(module.cmp0)
ftrxs = filterIncluTrxs(ftrxs)

ntrxs = []
ftrxsD = module.clusterOverTrxs(ftrxs)
for loci in ftrxsD.keys():
	ntrxs += module.renameTrxs2(loci, ftrxsD[loci])

totalTrxNum = len(trxs); ntotalTrxNum = len(ntrxs)

outputFile = open(outputdir + '/logs/transcripts_' + chr + '.filter.logs', 'w')
outputFile.write(chr + '\t' + str(totalTrxNum) + '\t' + str(ntotalTrxNum))
outputFile.close()

newGtf = dict(); newGtf[chr] = ntrxs
outputGtf = outputdir + '/transcripts_' + chr + '.tss.cps.filter.gtf'
module.writeGtf(newGtf, outputGtf)
