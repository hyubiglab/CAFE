#!/share/apps/python/bin/python
import sys, os

import config as conf
import data as data
import module as module

type = conf.tss_type

assembly = sys.argv[1]
gtfFile = sys.argv[2]
chr = sys.argv[3]
outputdir = sys.argv[4]
#type = sys.argv[5]

if assembly == 'hg19':
	fpseqAnno = data.hm_fpseqAll
	fpseqIntr = conf.hm_fpseqInterval; fpseqThr = conf.hm_fpseqThreshold
elif assembly == 'mm9':
	fpseq = data.ms_fpseqAll
	fpseqIntr = conf.ms_fpseqInterval; fpseqThr = conf.ms_fpseqThreshold

#if type == 'all':
#	inputdir = '/'.join(gtfFile.split('/')[:-1]) + '/pos/'
#	fpseqAnno = os.listdir(inputdir)
#	fpseqAnno = filter(lambda x: chr + '.tss.pos' in x, fpseqAnno)[0]
#	fpseqAnno = inputdir + fpseqAnno

def flankingIntr(geneid, firstExon, trxToCheck, interval, sense):
	trxToCheck = filter(lambda x: x.sense() == sense, trxToCheck)
	if sense == '+':
		trxToCheck = filter(lambda x: firstExon.start() - x.end() > 0 and (firstExon.start() - x.end()) + 50 < interval, trxToCheck) ##
		trxToCheck = filter(lambda x: not len(x.exons()) == 1, trxToCheck)
#		trxToCheck = filter(lambda x: not geneid == x.geneid(), trxToCheck) ##
		trxToCheck.reverse()
		if len(trxToCheck) > 0: interval = min(map(lambda x: (firstExon.start() - x.end()) + 50, trxToCheck)) ##
	else: #elif sense == '-':
		trxToCheck = filter(lambda x: x.start() - firstExon.end() > 0 and (x.start() - firstExon.end()) + 50 < interval, trxToCheck) ##
		trxToCheck = filter(lambda x: not len(x.exons()) == 1, trxToCheck)
#		trxToCheck = filter(lambda x: not geneid == x.geneid(), trxToCheck) ##
		if len(trxToCheck) > 0: interval = min(map(lambda x: (x.start() - firstExon.end()) + 50, trxToCheck)) ##
	return interval

def newTrx(trx, x):
	chr = trx.chr(); sense = trx.sense(); exons = trx.exons()
	nexons = []; nexon = ''

	if sense == '+': #forward strand
		start = x; end = trx.end()
		if len(exons) == 1: #single-exon transcript
			nexon = module.exon(chr, x, exons[-1].end(), sense)
			nexons.append(nexon)
		else: #multi-exon transcript
			for i in xrange(1, len(exons)):
				if exons[-i-1].end() > x:
					nexons.append(exons[-i])
					if i == len(exons)-1:
						nexon = module.exon(chr, x, exons[0].end(), sense)
						nexons.append(nexon)
				else:
					nexon = module.exon(chr, x, exons[-i].end(), sense)
					nexons.append(nexon)
					break
	else: #reverse strand
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

	nexons.sort(module.cmp0)
#	ntrxid = trx.trxid() + '_TSS.' + str(x)
	ntrxid = trx.trxid().split('_TSS')[0] + '_TSS.' + str(x)
	ntrx = module.transcript(chr, trx.geneid(), ntrxid, nexons[0].start(), nexons[-1].end(), nexons, trx.sense())
	ntrx.setFpPos(x)
	return ntrx

def filterOverTrxs(trxs):
	ntrxs = []; trxNum = len(trxs)
	for i in xrange(trxNum):
		trx = trxs[i]; flag = True
#		if trx.getFpPos() < 0:
		if trx.getFpPos() < 0 or len(trx.exons()) == 1:
			trxToCheck = trxs[max(0, i-20):i]
			trxToCheck += trxs[i+1:min(trxNum, i+20)]
			for trxT in trxToCheck:
				if trxT.getFpPos() > 0 and module.overlappedTrxs(trx, trxT):
					flag = False; break
			del trxToCheck
			if trx.getFpPos() < 0 and len(trx.exons()) == 1: flag = False
#			if trx.trxid().find('LM') > 0: flag = True
			if flag: ntrxs += [trx]
		else:
			if len(trx.exons()) == 1:
				trxToCheck = trxs[max(0, i-20):i]
				trxToCheck += trxs[i+1:min(trxNum, i+20)]
				for trxT in trxToCheck:
					if trxT.getFpPos() == trx.getFpPos() and trxT.exonNum() > 1:
						flag = False; break
				del trxToCheck
				if flag: ntrxs += [trx]
			else: ntrxs += [trx]
	return ntrxs

def updating_tss(trxs, fpseqAnno, fpseqIntr, fpseqThr, type):
	cSites = dict(); cSites['+'] = dict(); cSites['-'] = dict()
	uSites = dict(); uSites['+'] = dict(); uSites['-'] = dict()
#	if type == 'all':
#		cSiteD = dict(); cSiteD[chr] = []
#		fpseqAnno = open(fpseqAnno)
#		for line in fpseqAnno:
#			line = line.split('\t')
#			if len(line) == 4: cSiteD[chr].append([int(line[1]), float(line[2]), line[3].strip()])
#	else: cSiteD = module.readingAnno(fpseqAnno, 'cage')
	cSiteD = module.readingAnno(fpseqAnno, 'cage')
	for cSite in cSiteD[chr]: cSites[cSite[2]][cSite[0]] = ''

	totalTrxNum = 0; correctedTrxNum = 0; ntotalTrxNum = 0
	ntrxs = []; trxNum = len(trxs); trxs.sort(module.cmp0)
	for i in xrange(trxNum):
		trx = trxs[i]
		exons = trx.exons(); sense = trx.sense()

		interval = fpseqIntr #upstream interval
		ntrx = []; trxToCheck = []; fpEnds = []; exonToCheck = exons
		if sense == '+':
			if len(exonToCheck) == 0: exonToCheck.append(module.exon(exons[0].chr(), exons[0].start()-1, exons[0].start()-1, exons[0].sense()))
			firstExon = exonToCheck[0]; lastExon = exonToCheck[-1]
			trxToCheck = filter(lambda x: x.sense() == sense, trxs[:i])
		else: #elif sense == '-':
			if len(exonToCheck) == 0: exonToCheck.append(module.exon(exons[-1].chr(), exons[-1].end()+1, exons[-1].end()+1, exons[-1].sense()))
			firstExon = exonToCheck[-1]; lastExon = exonToCheck[0]
			trxToCheck = filter(lambda x: x.sense() == sense, trxs[i+1:])
		if len(trxToCheck) > 0: interval = flankingIntr(trx.geneid(), firstExon, trxToCheck, interval, sense)

		fpEnds, ufpEnds = module.getTSS(firstExon, lastExon, exonToCheck, cSiteD, interval, fpseqThr, type)
		for fpEnd in fpEnds: uSites[firstExon.sense()][fpEnd] = ''
		for ufpEnd in ufpEnds: uSites[firstExon.sense()][ufpEnd] = ''
		
		if len(fpEnds) > 0:
			ntrx = map(lambda x: newTrx(trx, x), fpEnds)
			fpEnds.sort()
			j = -1
			if sense == '+' and fpEnds[0] > firstExon.end(): j = 1 ##
			elif sense == '-' and fpEnds[-1] < firstExon.start(): j = 1 ##
			if j > 0:
				ntrxid = trx.trxid() + 'LM'
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
	totalTssNum = sum(map(lambda x: len(cSites[x].keys()), cSites.keys()))
	assignTssNum = sum(map(lambda x: len(uSites[x].keys()), uSites.keys()))

	# printing output
#	outputFile = open(outputdir + '/logs/transcripts_' + chr + '.tss.' + type + '.logs', 'w')
	outputFile = open(outputdir + '/logs/transcripts_' + chr + '.tss.logs', 'w')
	outputFile.write(chr + '\t' + str(totalTrxNum) + '\t' + str(correctedTrxNum) + '\t' + str(ntotalTrxNum) + '\t' + str(totalTssNum) + '\t' + str(assignTssNum))
	outputFile.close()

	newGtf = dict(); newGtf[chr] = ntrxs
#	outputGtf = outputdir + '/transcripts_' + chr + '.tss.' + type + '.gtf'
	outputGtf = outputdir + '/transcripts_' + chr + '.tss.gtf'
	module.writeGtf(newGtf, outputGtf)

#	if not type == 'all':
#		nSiteD = dict(); nSiteD[chr] = []
#		for cSite in cSiteD[chr]:
#			if not cSite[0] in uSites[cSite[2]].keys(): nSiteD[chr].append(cSite)
#		if not os.path.exists(outputdir + '/pos/'): os.makedirs(outputdir + '/pos/')
#		outputFile2 = open(outputdir + '/pos/transcripts_' + chr + '.tss.pos', 'w')
#		for nSite in nSiteD[chr]:
#			outputFile2.write(chr + '\t' + str(nSite[0]) + '\t' + str(nSite[1]) + '\t' + str(nSite[2]) + '\n')
#		outputFile2.close()

gtfFile = module.getGtf(gtfFile)
trxs = gtfFile[chr]
updating_tss(trxs, fpseqAnno, fpseqIntr, fpseqThr, type)
