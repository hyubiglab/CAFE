#!/share/apps/python/bin/python
import sys, os

import config as conf
import data as data
import module as module

type = conf.cps_type

assembly = sys.argv[1]
gtfFile = sys.argv[2]
chr = sys.argv[3]
outputdir = sys.argv[4]
#type = sys.argv[5]

if assembly == 'hg19':
	tpseqAnno = data.hm_tpseqAll
	tpseqIntr = conf.hm_tpseqInterval; tpseqThr = conf.hm_tpseqThreshold
elif assembly == 'mm9':
	tpseq = data.ms_tpseqAll
	tpseqIntr = conf.ms_tpseqInterval; tpseqThr = conf.ms_tpseqThreshold

#if type == 'all':
#	inputdir = '/'.join(gtfFile.split('/')[:-1]) + '/pos/'
#	tpseqAnno = os.listdir(inputdir)
#	tpseqAnno = filter(lambda x: chr + '.cps.pos' in x, tpseqAnno)[0]
#	tpseqAnno = inputdir + tpseqAnno

def flankingIntr(geneid, lastExon, trxToCheck, interval, sense):
	trxToCheck = filter(lambda x: x.sense() == sense, trxToCheck)
	if sense == '+':
		trxToCheck = filter(lambda x: x.start() - lastExon.end() > 0 and (x.start() - lastExon.end()) + 100 < interval, trxToCheck) ##
		trxToCheck = filter(lambda x: not len(x.exons()) == 1, trxToCheck)
#		trxToCheck = filter(lambda x: not geneid == x.geneid(), trxToCheck) ##
		if len(trxToCheck) > 0: interval = min(map(lambda x: (x.start() - lastExon.end()) + 100, trxToCheck)) ##
	else: #elif sense == '-':
		trxToCheck = filter(lambda x: lastExon.start() - x.end() > 0 and (lastExon.start() - x.end()) + 100 < interval, trxToCheck) ##
		trxToCheck = filter(lambda x: not len(x.exons()) == 1, trxToCheck)
#		trxToCheck = filter(lambda x: not geneid == x.geneid(), trxToCheck) ##
		trxToCheck.reverse()
		if len(trxToCheck) > 0: interval = min(map(lambda x: (lastExon.start() - x.end()) + 100, trxToCheck)) ##
	return interval

def newTrx(trx, x):
	chr = trx.chr(); sense = trx.sense(); exons = trx.exons()
	nexons = []; nexon = ''

	if sense == '+': #forward strand
		start = trx.start(); end = x
		if len(exons) == 1: #single-exon transcript
			nexon = module.exon(chr, exons[-1].start(), x, sense)
			nexons.append(nexon)
		else: #multi-exon transcript
			for i in xrange(1, len(exons)):
				if exons[i].start() < x:
					nexons.append(exons[i-1])
					if i == len(exons) - 1:
						nexon = module.exon(chr, exons[-1].start(), x, sense)
						nexons.append(nexon)
				else:
					nexon = module.exon(chr, exons[i-1].start(), x, sense)
					nexons.append(nexon)
					break
	else: #reverse strand
		start = x; end = trx.end()
		if len(exons) == 1: #single-exon transcript
			nexon = module.exon(chr, x, exons[-1].end(), sense)
			nexons.append(nexon)
		else: #multi-exon transcript
			for i in xrange(1, len(exons)):
				if exons[-i-1].end() > x:
					nexons.append(exons[-i])
					if i == len(exons) - 1:
						nexon = module.exon(chr, x, exons[0].end(), sense)
						nexons.append(nexon)
				else:
					nexon = module.exon(chr, x, exons[-i].end(), sense)
					nexons.append(nexon)
					break

	nexons.sort(module.cmp0)
#	ntrxid = trx.trxid() + '_CPS.' + str(x)
	ntrxid = trx.trxid().split('_CPS')[0] + '_CPS.' + str(x)
	ntrx = module.transcript(chr, trx.geneid(), ntrxid, nexons[0].start(), nexons[-1].end(), nexons, trx.sense())
	ntrx.setTpPos(x)
	return ntrx

def filterOverTrxs(trxs):
	ntrxs = []; trxNum = len(trxs)
	for i in xrange(trxNum):
		trx = trxs[i]; flag = True
#		if trx.getTpPos() < 0:
		if trx.getTpPos() < 0 or len(trx.exons()) == 1:
			trxToCheck = trxs[max(0, i-20):i]
			trxToCheck += trxs[i+1:min(trxNum, i+20)]
			for trxT in trxToCheck:
				if trxT.getTpPos() > 0 and module.overlappedTrxs(trx, trxT):
					flag = False; break
			del trxToCheck
			if trx.getTpPos() < 0 and len(trx.exons()) == 1: flag = False
#			if trx.trxid().find('RM') > 0: flag = True
			if flag: ntrxs += [trx]
		else:
			if len(trx.exons()) == 1:
				trxToCheck = trxs[max(0, i-20):i]
				trxToCheck += trxs[i+1:min(trxNum, i+20)]
				for trxT in trxToCheck:
					if trxT.getTpPos() == trx.getTpPos() and trxT.exonNum() > 1:
						flag = False; break
				del trxToCheck
				if flag: ntrxs += [trx]
			else: ntrxs += [trx]
	return ntrxs

def updating_cps(trxs, tpseqAnno, tpseqIntr, tpseqThr, type):
	cSites = dict(); cSites['+'] = dict(); cSites['-'] = dict()
	uSites = dict(); uSites['+'] = dict(); uSites['-'] = dict()
#	if type == 'all':
#		cSiteD = dict(); cSiteD[chr] = []
#		tpseqAnno = open(tpseqAnno)
#		for line in tpseqAnno:
#			line = line.split('\t')
#			if len(line) == 4: cSiteD[chr].append([int(line[1]), float(line[2]), line[3].strip()])
#	else: cSiteD = module.readingAnno(tpseqAnno, 'polya')
	cSiteD = module.readingAnno(tpseqAnno, 'polya')
	for cSite in cSiteD[chr]: cSites[cSite[2]][cSite[0]] = ''

	totalTrxNum = 0; correctedTrxNum = 0; ntotalTrxNum = 0
	ntrxs = []; trxNum = len(trxs); trxs.sort(module.cmp0)
	for i in xrange(trxNum):
		trx = trxs[i]
		exons = trx.exons(); sense = trx.sense()

		interval = tpseqIntr #downstream interval
		ntrx = []; trxToCheck = []; tpEnds = []; exonToCheck = exons
		if sense == '+':
			if len(exonToCheck) == 0: exonToCheck.append(module.exon(exons[-1].chr(), exons[-1].end()+1, exons[-1].end()+1, exons[-1].sense()))
			firstExon = exonToCheck[0]; lastExon = exonToCheck[-1]
			trxToCheck = filter(lambda x: x.sense() == sense, trxs[i+1:])
		else: #elif sense == '-':
			if len(exonToCheck) == 0: exonToCheck.append(module.exon(exons[0].chr(), exons[0].end()-1, exons[0].end()-1, exons[0].sense()))
			firstExon = exonToCheck[-1]; lastExon = exonToCheck[0]
			trxToCheck = filter(lambda x: x.sense() == sense, trxs[:i])
		if len(trxToCheck) > 0: interval = flankingIntr(trx.geneid(), lastExon, trxToCheck, interval, sense)

		tpEnds, utpEnds = module.getCPS(firstExon, lastExon, exonToCheck, cSiteD, interval, tpseqThr, type)
		for tpEnd in tpEnds: uSites[lastExon.sense()][tpEnd] = ''
		for utpEnd in utpEnds: uSites[lastExon.sense()][utpEnd] = ''

		if len(tpEnds) > 0:
			ntrx = map(lambda x: newTrx(trx, x), tpEnds)
			tpEnds.sort()
			j = -1
			if sense == '+' and tpEnds[-1] < lastExon.start(): j = 1 ##
			elif sense == '-' and tpEnds[0] > lastExon.end(): j = 1 ##
			if j > 0:
				ntrxid = trx.trxid() + '_RM'
				ntrxs.append(module.transcript(trx.chr(), trx.geneid(), ntrxid, trx.start(), trx.end(), exons, sense))

		if len(ntrx) > 0: ntrxs += ntrx; correctedTrxNum += 1
		else: ntrxs += [trx]

	# filter transcripts
	ntrxs.sort(module.cmp0)
	ntrxs = module.filterSameTrxs(ntrxs)
	ntrxs = module.filterNoneTrxs(ntrxs)
	ntrxs = module.checkProperTrxs(ntrxs)
#	if type == 'all': ntrxs = filterOverTrxs(ntrxs)
	ntrxs = filterOverTrxs(ntrxs)

	totalTrxNum += len(trxs); ntotalTrxNum += len(ntrxs)
	totalCpsNum = sum(map(lambda x: len(cSites[x].keys()), cSites.keys()))
	assignCpsNum = sum(map(lambda x: len(uSites[x].keys()), uSites.keys()))

	# printing output
#	outputFile = open(outputdir + '/logs/transcripts_' + chr + '.tss.cps' + type + '.logs', 'w')
	outputFile = open(outputdir + '/logs/transcripts_' + chr + '.cps.logs', 'w')
	outputFile.write(chr + '\t' + str(totalTrxNum) + '\t' + str(correctedTrxNum) + '\t' + str(ntotalTrxNum) + '\t' + str(totalCpsNum) + '\t' + str(assignCpsNum))
	outputFile.close()

	newGtf = dict(); newGtf[chr] = ntrxs
#	outputGtf = outputdir + '/transcripts_' + chr + '.tss.cps.' + type + '.gtf'	
	outputGtf = outputdir + '/transcripts_' + chr + '.tss.cps.gtf'
	module.writeGtf(newGtf, outputGtf)
	
#	if not type == 'all':
#		nSiteD = dict(); nSiteD[chr] = []
#		for cSite in cSiteD[chr]:
#			if not cSite[0] in uSites[cSite[2]].keys(): nSiteD[chr].append(cSite)
#		if not os.path.exists(outputdir + '/pos/'): os.makedirs(outputdir + '/pos/')
#		outputFile2 = open(outputdir + '/pos/transcripts_' + chr + '.cps.pos', 'w')
#		for nSite in nSiteD[chr]:
#			outputFile2.write(chr + '\t' + str(nSite[0]) + '\t' + str(nSite[1]) + '\t' + str(nSite[2]) + '\n')
#		outputFile2.close()

gtfFile = module.getGtf(gtfFile)
trxs = gtfFile[chr]
updating_cps(trxs, tpseqAnno, tpseqIntr, tpseqThr, type)
