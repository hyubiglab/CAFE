#!/share/apps/python/bin/python
import sys, os, copy
from subprocess import Popen, PIPE, STDOUT

import config as conf
import data as data
import module as module
import sam as sam

type = conf.exj_type
cutoff = float(conf.exj_cutoff)

assembly = sys.argv[1]
gtfFile = sys.argv[2]
rpds = sys.argv[3]
stranded = sys.argv[4]
chr = sys.argv[5]
outputdir = sys.argv[6]

if assembly == 'hg19':
	fasta = data.hm_fasta
	intronmin = conf.hm_intronmin; intronmax = conf.hm_intronmax
	fpsite_dic = module.readingAnno(data.hm_fpseqAll, 'cage')
	tpsite_dic = module.readingAnno(data.hm_tpseqAll, 'polya')
elif assembly == 'mm9':
	fasta = data.ms_fasta
	intronmin = conf.ms_intronmin; intronmax = conf.ms_intronmax
	fpsite_dic = module.readingAnno(data.ms_fpseqAll, 'cage')
	tpsite_dic = module.readingAnno(data.ms_tpseqAll, 'polya')

seq = ''
if type == 'signal':
	seq_dic = module.getSeq_dictbyChr(fasta, chr); seq = seq_dic[chr]
#	score5 = data.datadir[:-5] + '/build/MaxEntScan/score5.pl'
#	score3 = data.datadir[:-5] + '/build/MaxEntScan/score3.pl'
	score5 = data.datadir[:-5] + '/build/MaxEntScan/score5.pl ' + data.datadir[:-5] + '/build/MaxEntScan/'
	score3 = data.datadir[:-5] + '/build/MaxEntScan/score3.pl ' + data.datadir[:-5] + '/build/MaxEntScan/'

noutputdir = outputdir + '/bam/'
if not os.path.exists(noutputdir): os.makedirs(noutputdir)
merged_bamfile = noutputdir + '/' + stranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + '.merged.bam'
sam.bamCat([rpds, stranded], merged_bamfile)

def checkIntrons(intron, i, sortedstarts, trxs):
	if i >= len(sortedstarts)-10: trxToCheck = sortedstarts[i-10:len(sortedstarts)]
	elif i < 10: trxToCheck = sortedstarts[0:i+10]
	else: trxToCheck = sortedstarts[i-10:i+10]
	for geneid in trxToCheck:
		for isoTrx in trxs[geneid]:
			for intron2 in isoTrx.introns():
				if intron.overlaps(intron2):
					return True
	return False

def checkSites(firstExon, lastExon, cSiteD):
	cSites = []
	for cSite in cSiteD[firstExon.chr()]:
		if (int(cSite[0]) >= firstExon.start()) and (int(cSite[0]) <= lastExon.end()) and (cSite[2] == firstExon.sense()):
			cSites.append(cSite)
			break
	if len(cSites) == 0: return True
	else: return False

#def findSignal(seq, signal, indexList, offset): #recursive function
#	if seq.find(signal) == -1:
#		return indexList
#	else:
#		offset = seq.index(signal)+offset
#		indexList.append(offset)
#		seq = seq[seq.index(signal)+1:]
#		return findSignal(seq, signal, indexList, offset+1)

def findSignal(seq, signal, indexList, offset): #loop function
	for i in xrange(seq.count(signal)):
		offset = seq.index(signal)+offset
		indexList.append(offset)
		seq = seq[seq.index(signal)+1:]
		offset += 1
	return indexList

def MaxEntScan(seq, type):
	if type == '5ss': commandBase = 'perl ' + score5 + ' ' + str(seq)
	elif type == '3ss': commandBase = 'perl ' + score3 + ' ' + str(seq)
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	score = p.stdout.read(); score = score.strip()
	return score

def updating_exj(trxs, sense, seq, type = 'reads'):
	genestarts = dict()
	# sort the gene dictionary by start position
	for geneid in trxs.keys(): genestarts[geneid] = trxs[geneid][0].start()
	starts = [(v, k) for k, v in genestarts.items()]
	sortedstarts = [k for v, k in sorted(starts)]
	totalJunctions = 0; readJunctions = 0; signalJunctions = 0
	# reconstruct the transfrags based on exon-junctions
	ntrxs = copy.deepcopy(trxs[sortedstarts[0]])
	for i in xrange(1, len(sortedstarts)):
		geneid = sortedstarts[i]

		def finding_exj(ntrxs, geneid, sortedstarts, trxs, totalJunctions, readJunctions, signalJunctions):
			index_dic = dict(); nntrxs = []
			for i in xrange(0, len(ntrxs)):
				index_dic[i] = []; j = 0
				newTrx = ntrxs[i]
				for isoTrx in trxs[geneid]:
					if not newTrx.geneid() == isoTrx.geneid():
						exon1 = newTrx.exons()[-1]; exon2 = isoTrx.exons()[0]
#						distance = exon2.start() - exon1.end() + 1 ##
						distance = exon2.end() - exon1.start() + 1 ##
						if intronmin <= distance <= intronmax:
							intron = module.exon(chr, exon1.end(), exon2.start(), exon1.sense()); totalJunctions += 1
							if not checkIntrons(intron, sortedstarts.index(geneid), sortedstarts, trxs):
								flag, nexons = module.updateIntrons(merged_bamfile, exon1, exon2)
								if flag:
#									ngeneid = newTrx.geneid() + '_' + isoTrx.geneid().split('.')[1] ##
#									ntrxid = newTrx.trxid() + '_' + '.'.join(isoTrx.trxid().split('.')[1:]) ##
									ngeneid = newTrx.geneid() + '_' + isoTrx.geneid(); ntrxid = newTrx.trxid() + '_' + isoTrx.trxid() ##
									if len(newTrx.exons()) == 1 and len(isoTrx.exons()) == 1: nexons1 = nexons; nexons2 = []
									elif len(newTrx.exons()) > 1 and len(isoTrx.exons()) == 1: nexons1 = newTrx.exons()[:-1]; nexons2 = nexons
									elif len(newTrx.exons()) == 1 and len(isoTrx.exons()) > 1: nexons1 = nexons; nexons2 = isoTrx.exons()[1:]
									else: nexons1 = newTrx.exons()[:-1]; nexons1 += nexons; nexons2 = isoTrx.exons()[1:]
									index_dic[i].append(copy.deepcopy(newTrx))
									index_dic[i][-1].update(ngeneid, ntrxid, newTrx.start(), isoTrx.end(), nexons1, nexons2)
									isoTrx.setExj(True); j += 1; readJunctions += 1
								else:
									if len(newTrx.exons()) == 1 or len(isoTrx.exons()) == 1:
										if checkSites(exon1, exon2, fpsite_dic) and checkSites(exon1, exon2, tpsite_dic):
											flag = True
									else:
										if checkSites(exon1, exon2, fpsite_dic) or checkSites(exon1, exon2, tpsite_dic):
											flag = True
									if flag:
										# find potential exon-junctions
#										exon1_seq = seq[exon1.start()+1:exon1.end()+10]; exon2_seq = seq[exon2.start()-11:exon2.end()-1] ##
										exon1_seq = seq[exon1.start()+1:exon1.end()+100]; exon2_seq = seq[exon2.start()-101:exon2.end()-1] ##
										exon1_dict = dict(); exon2_dict = dict()
										if sense == '+':
											donor = 'GT'; acceptor = 'AG'
											donor_index = findSignal(exon1_seq, donor, [], exon1.start()+1)
											for index in donor_index:
												exon1_dict[index] = seq[index-3:index+6]
#											acceptor_index = findSignal(exon2_seq, acceptor, [], exon2.start()-11) ##
											acceptor_index = findSignal(exon2_seq, acceptor, [], exon2.start()-101) ##
											for index in acceptor_index:
												exon2_dict[index] = seq[index-18:index+5]
										elif sense == '-':
											donor = 'AC'; acceptor = 'CT'
#											donor_index = findSignal(exon2_seq, donor, [], exon2.start()-11) ##
											donor_index = findSignal(exon2_seq, donor, [], exon2.start()-101) ##
											for index in donor_index:
												exon1_dict[index] = module.reverseComp(seq[index-4:index+5])
											acceptor_index = findSignal(exon1_seq, acceptor, [], exon1.start()+1)
											for index in acceptor_index:
												exon2_dict[index] = module.reverseComp(seq[index-3:index+20])
										if len(exon1_dict.keys()) > 0 and len(exon2_dict.keys()) > 0:
											for fss_index in exon1_dict.keys():
												if exon1_dict[fss_index].find('N') < 0:
													exon1_dict[fss_index] = float(MaxEntScan(exon1_dict[fss_index], '5ss'))
												else: exon1_dict[fss_index] = float(-20.0)
											for tss_index in exon2_dict.keys():
												if exon2_dict[tss_index].find('N') < 0:
													exon2_dict[tss_index] = float(MaxEntScan(exon2_dict[tss_index], '3ss'))
												else: exon2_dict[tss_index] = float(-20.0)
											fss_max = max(exon1_dict, key=exon1_dict.get); tss_max = max(exon2_dict, key=exon2_dict.get)
											if exon1_dict[fss_max] >= cutoff and exon2_dict[tss_max] >= cutoff:
#												ngeneid = newTrx.geneid() + '_' + isoTrx.geneid().split('.')[1] ##
#												ntrxid = newTrx.trxid() + '_' + '.'.join(isoTrx.trxid().split('.')[1:]) ##
												ngeneid = newTrx.geneid() + '_' + isoTrx.geneid(); ntrxid = newTrx.trxid() + '_' + isoTrx.trxid() ##
												nnewTrx = copy.deepcopy(newTrx); nisoTrx = copy.deepcopy(isoTrx)
												if sense == '+': nnewTrx.exons()[-1:][0].newEnd(fss_max); nisoTrx.exons()[0].newStart(tss_max+3)
												elif sense == '-': nnewTrx.exons()[-1:][0].newEnd(tss_max); nisoTrx.exons()[0].newStart(fss_max+3)
												index_dic[i].append(copy.deepcopy(nnewTrx)); index_dic[i][-1].update(ngeneid, ntrxid, nnewTrx.start(), nisoTrx.end(), nnewTrx.exons(), nisoTrx.exons())
												isoTrx.setExj(True); j += 1; signalJunctions += 1
				if j == 0: index_dic[i].append(copy.deepcopy(newTrx))
			index_list = index_dic.keys(); index_list.sort()
			for x in index_list:
				for y in index_dic[x]: nntrxs.append(y)
			for z in trxs[geneid]:
				if not z.getExj(): nntrxs.append(z)
			return nntrxs, totalJunctions, readJunctions, signalJunctions

		ntrxs, totalJunctions, readJunctions, signalJunctions = finding_exj(ntrxs, geneid, sortedstarts, trxs, totalJunctions, readJunctions, signalJunctions)
	return ntrxs, totalJunctions, readJunctions, signalJunctions

gtfFile = module.getGtf(gtfFile)
trxs = gtfFile[chr]

ftrxs = dict(); rtrxs = dict()
trxNum = len(trxs)
for i in xrange(trxNum):
	trx = trxs[i]
	if trx.sense() == '+':
		if not ftrxs.has_key(trx.geneid()):
			ftrxs[trx.geneid()] = []
		ftrxs[trx.geneid()].append(trx)
	else: #elif trx.sense() == '-':
		if not rtrxs.has_key(trx.geneid()):
			rtrxs[trx.geneid()] = []
		rtrxs[trx.geneid()].append(trx)

# forward strand
nftrxs, ftotalJunctions, freadJunctions, fsignalJunctions = updating_exj(ftrxs, '+', seq, type)

# reverse strand
nrtrxs, rtotalJunctions, rreadJunctions, rsignalJunctions = updating_exj(rtrxs, '-', seq, type)

ntrxs = []; ntrxs += nftrxs; ntrxs += nrtrxs
ntrxs = module.filterSameTrxs(ntrxs)

ntrxs.sort(module.cmp0)
ntrxs = module.filterIncluTrxs(ntrxs)

totalTrxNum = len(ntrxs); totalJunctions = ftotalJunctions + rtotalJunctions
readJunctions = freadJunctions + rreadJunctions; signalJunctions = fsignalJunctions + rsignalJunctions
outputFile = open(outputdir + '/logs/transcripts_' + chr + '.exj.logs', 'w')
outputFile.write(chr + '\t' + str(trxNum) + '\t' + str(totalTrxNum) + '\t' + str(totalJunctions) + '\t' + str(readJunctions) + '\t' + str(signalJunctions))
outputFile.close()

newGtf = dict(); newGtf[chr] = ntrxs
outputGtf = outputdir + '/transcripts_' + chr + '.exj.gtf'
module.writeGtf(newGtf, outputGtf)
