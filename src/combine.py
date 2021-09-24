#!/share/apps/python/bin/python
import sys, os
sys.setrecursionlimit(10000)

import module as module

assembly = sys.argv[1]
gtfFiles = sys.argv[2]
chr = sys.argv[3]
outputdir = sys.argv[4]

def sameLocus(locus, otherLocus):
	if locus.chr() == otherLocus.chr() and locus.sense() == otherLocus.sense() and locus.start() == otherLocus.start() and locus.end() == otherLocus.end(): return True
	else: return False

def overLocus(locus, otherLocus):
	if locus.chr() != otherLocus.chr(): return False
	elif not(locus.sense() == '.' or \
		otherLocus.sense() == '.' or \
		locus.sense() == otherLocus.sense()): return False
	elif locus.start() > otherLocus.end() or otherLocus.start() > locus.end(): return False
	else: return True

def updateOverlappedTrx(exons, exons2, introns, introns2):
	introns += introns2
	introns_start = map(lambda x: x.start(), introns)
	introns_end = map(lambda y: y.end(), introns)
	introns_start = dict.fromkeys(introns_start).keys(); introns_start.sort()
	introns_end = dict.fromkeys(introns_end).keys(); introns_end.sort()
	if not len(introns_start) == len(introns_end):
		return False, []
	else:
		nintrons = []
		for n in xrange(len(introns_start)):
			nintrons.append([introns_start[n], introns_end[n]])
		nexons = []
		for i in xrange(len(nintrons)):
			if i == 0: nexon = module.exon(exons[0].chr(), min(exons[0].start(), exons2[0].start()), nintrons[i][0]-1, exons[0].sense())
			else: nexon = module.exon(exons[0].chr(), nintrons[i-1][1]+1, nintrons[i][0]-1, exons[0].sense())
			nexons.append(nexon)
			if i == len(nintrons) - 1:
				nexon = module.exon(exons[0].chr(), nintrons[i][1]+1, max(exons[-1].end(), exons2[-1].end()), exons[0].sense())
		return True, nexons

def combineOverTrxs(loci, trxs):
	trxs = module.renameTrxs(loci, trxs)
	if len(trxs) == 1: return trxs
	else:
		trxs.sort(module.cmp0)
		ntrxs = []; flag = False
		for x in xrange(len(trxs)):
			for y in xrange(len(trxs)):
				if x < y:
					nexons = []; flag2 = False
					trx = trxs[x]; trx2 = trxs[y]
					exons = trx.exons(); exons2 = trx2.exons()

					if len(exons) == 1 and len(exons2) == 1:
						if exons[0].overlaps(exons2[0]):
							nexons.append(module.exon(trx.chr(), min(exons[0].start(), exons2[0].start()), max(exons[-1].end(), exons2[-1].end()), trx.sense()))
							flag = True; flag2 = True
					elif len(exons) == 1:
						nexons += exons2
						if exons[0].overlaps(exons2[0]):
							nexons[0] = module.exon(trx.chr(), min(exons[0].start(), exons2[0].start()), exons2[0].end(), trx.sense())
							flag = True; flag2 = True
						if exons[0].overlaps(exons2[-1]):
							nexons[-1] = module.exon(trx.chr(), exons2[-1].start(), max(exons[-1].end(), exons2[-1].end()), trx.sense())
							flag = True; flag2 = True
					elif len(exons2) == 1:
						nexons += exons
						if exons2[0].overlaps(exons[0]):
							nexons[0] = module.exon(trx.chr(), min(exons[0].start(), exons2[0].start()), exons[0].end(), trx.sense())
							flag = True; flag2 = True
						if exons2[0].overlaps(exons[-1]):
							nexons[-1] = module.exon(trx.chr(), exons[-1].start(), max(exons[-1].end(), exons2[-1].end()), trx.sense())
							flag = True; flag2 = True
					else:
						exonNum = len(exons); exonNum2 = len(exons2)
						introns = trx.introns(); introns2 = trx2.introns()
						intronNum = len(introns); intronNum2 = len(introns2)

						sameIntron = 0; overIntron = 0; noverIntron = 0
						sameIntron2 = 0; overIntron2 = 0; noverIntron2 = 0
						for intron in introns:
							iflag = False
							for intron2 in introns2:
								if sameLocus(intron, intron2):
									sameIntron += 1; iflag = True
									break
								elif overLocus(intron, intron2):
									overIntron += 1; iflag = True
									break
							if not iflag: noverIntron += 1
						for intron2 in introns2:
                                                        iflag2 = False
                                                        for intron in introns:
                                                                if sameLocus(intron2, intron):
                                                                        sameIntron2 += 1; iflag2 = True
                                                                        break
                                                                elif overLocus(intron2, intron):
                                                                        overIntron2 += 1; iflag2 = True
                                                                        break
                                                        if not iflag2: noverIntron2 += 1

						if overIntron == 0 and overIntron2 == 0:
							if intronNum == sameIntron and intronNum2 == sameIntron2:
								nexons += exons
								nexons[0] = module.exon(trx.chr(), min(exons[0].start(), exons2[0].start()), exons[0].end(), trx.sense())
								nexons[-1] = module.exon(trx.chr(), exons[-1].start(), max(exons[-1].end(), exons2[-1].end()), trx.sense())
								flag = True; flag2 = True
							elif intronNum == noverIntron and intronNum2 == noverIntron2:
								if exons2[0].start() <= exons[-1].end() < exons2[1].start() and exons[-2].end() < exons2[0].start() <= exons[-1].end(): ##
									nexons += exons
									nexons[-1] = module.exon(trx.chr(), exons[-1].start(), exons2[0].end(), trx.sense())
									nexons += exons2[1:]
									flag = True; flag2 = True
							else:
								eflag = False
								for intron in introns:
									for exon2 in exons2:
										if overLocus(intron, exon2):
											eflag = True
											break
									if eflag: break
								eflag2 = False
								for intron2 in introns2:
									for exon in exons:
										if overLocus(intron2, exon):
											eflag2 = True
											break
									if eflag2: break
								if not eflag and not eflag2:
									uflag, nexons = updateOverlappedTrx(exons, exons2, introns, introns2)
									if uflag:
										flag = True; flag2 = True

					if flag2:
						nexons.sort(module.cmp0)
						ntrxs.append(module.transcript(trx.chr(), 'Gene' + trx.sense() + trx.geneid().split(trx.sense())[-1] + '_' + trx2.geneid().split(trx.sense())[-1], 'Trx' + trx.sense() + trx.trxid().split(trx.sense())[-1] + '_' + trx2.trxid().split(trx.sense())[-1], nexons[0].start(), nexons[-1].end(), nexons, trx.sense()))
						trx.setUse(True); trx2.setUse(True)

		for rtrx in trxs:
			if not rtrx.getUse(): ntrxs.append(rtrx)
		ntrxs = module.filterSameTrxs(ntrxs)
		ntrxs = module.renameTrxs(loci, ntrxs)
		if flag: return combineOverTrxs(loci, ntrxs)
		else: return ntrxs

gtfFiles = gtfFiles.split(',')
trxs = module.getGtf(gtfFiles[0])[chr]
for i in xrange(1, len(gtfFiles)):
	trxs += module.getGtf(gtfFiles[i])[chr]

totalTrxNum = 0; totalTrxNum += len(trxs)
trxs = module.fillIntrons(assembly, trxs)
trxs = module.filterSameTrxs(trxs)

trxs.sort(module.cmp0)
trxs = module.filterIncluTrxs(trxs)

ftrxs = []; rtrxs = []
trxNum = len(trxs)
for i in xrange(trxNum):
	trx = trxs[i]
	if trx.sense() == '+': ftrxs.append(trx)
	elif trx.sense() == '-': rtrxs.append(trx)
#	else: print trx.trxid(), ' ', trx.coord()

# forward strand
ftrxs.sort(module.cmp0)
ftrxsD = module.clusterOverTrxs(ftrxs)
for ftrxs in ftrxsD.keys():
	ftrxsD[ftrxs] = combineOverTrxs(ftrxs, ftrxsD[ftrxs])

# reverse strand
rtrxs.sort(module.cmp0)
rtrxsD = module.clusterOverTrxs(rtrxs)
for rtrxs in rtrxsD.keys():
	rtrxsD[rtrxs] = combineOverTrxs(rtrxs, rtrxsD[rtrxs])

ntrxs = []
for ftrxs in ftrxsD.keys():
	for nftrx in ftrxsD[ftrxs]: ntrxs.append(nftrx)
for rtrxs in rtrxsD.keys():
	for nrtrx in rtrxsD[rtrxs]: ntrxs.append(nrtrx)

ntrxs = module.filterSameTrxs(ntrxs)

ntrxs.sort(module.cmp0)
ntrxs = module.filterIncluTrxs(ntrxs)

combinedTrxNum = 0; combinedTrxNum += len(ntrxs)
outputFile = open(outputdir + '/logs/transcripts_' + chr + '.logs', 'w')
outputFile.write(chr + '\t' + str(totalTrxNum) + '\t' + str(combinedTrxNum))
outputFile.close()

newGtf = dict(); newGtf[chr] = ntrxs
outputGtf = outputdir + '/transcripts_' + chr + '.gtf'
module.writeGtf(newGtf, outputGtf)
