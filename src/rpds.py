#!/share/apps/python/bin/python
import sys, os
import copy, math

import config as conf
import data as data
import module as module
import sam as sam

type = conf.rpds_type
order = conf.rpds_order

assembly = sys.argv[1]
unstranded = sys.argv[2]
stranded = sys.argv[3]
chr = sys.argv[4]
outputdir = sys.argv[5]

# matrix
if assembly == 'hg19': matrix_dir = data.hm_matrix
elif assembly == 'mm9': matrix_dir = data.ms_matrix

matrix_dic = dict()
if type == 'hmm':
	for x in xrange(1, order+1):
		matrix_dic[x] = dict()
		matrix = open(matrix_dir + 'markov_' + str(x))
		for line in matrix:
			line = line.split(':')
			matrix_dic[x][line[0]] = dict()
			for y in xrange(1, 5):
				matrix_dic[x][line[0]][line[y].split("'")[1]] = line[y+1].split(',')[0].split('}')[0]

def rpds(unstranded, stranded, chr):
	finalMappingReads_dic = dict(); remainMappingReads_dic = dict(); npExonJunctionPairs_dic = dict()
	npMappingReads_dic, npExonJunctionReads_dic, npReadPos_dic = module.getBamDict(unstranded, chr)
	spReadPos_dic = module.getBamDictPos(stranded, chr, 'fr')
	
	npReadPos = sorted(npReadPos_dic.keys()); spReadPos = sorted(spReadPos_dic.keys()); index = 0
	for npPos in npReadPos: #unstranded read position
		forwardReads = 0; reverseReads = 0; forwardRead = 0; reverseRead = 0
		if len(spReadPos[max(0, index):]) > 0:
			posList, distance = module.getAdjacentPos2(npPos, spReadPos[max(0, index):])
			if int(distance) <= 50: ##
				for i in range(len(posList)):
					if len(spReadPos_dic[posList[i]]) > 1: #both strand
						forwardRead = len(spReadPos_dic[posList[i]]['+']); reverseRead = len(spReadPos_dic[posList[i]]['-'])
					else:
						if spReadPos_dic[posList[i]].keys()[0] == '+':
							forwardRead = len(spReadPos_dic[posList[i]]['+']); reverseRead = 0
						else: #spReadPos_dic[posList[i]].keys()[0] == '-':
							forwardRead = 0; reverseRead = len(spReadPos_dic[posList[i]]['-'])
						forwardReads += forwardRead; reverseReads += reverseRead

		index = spReadPos.index(posList[0]) - 10
		fProb = 1; rProb = 1
		if type == 'hmm':
			if len(spReadPos[max(0, index):]) > 0:
				apReads_dict = module.getAdjacentPos(npPos, spReadPos[max(0, index):])
				apReads_list = apReads_dict.keys(); apReads_list.sort()
				sense_list = []; nsense_list = []
				if int(apReads_list[0]) <= 50: ##
					snReads_dict = dict(); snReads_dict[0] = []
					for snRead in apReads_dict[apReads_list[0]]: snReads_dict[0].append([snRead])
					for x in xrange(1, order):
						if len(spReadPos[max(0, index):]) > x and int(snReads_dict.keys()[0]) == x-1:
							snReads_dict = module.getAdjacentPos4(snReads_dict, spReadPos[max(0, index):], x)

					for snReads_pos in snReads_dict.values()[0]:
						nsense_list2 = []
						for snRead_pos in snReads_pos: nsense_list2 += [spReadPos_dic[snRead_pos].keys()]
						nsense_list.append(nsense_list2)
					for m in nsense_list:
						nsense_list3 = []
						for n in xrange(0, len(m)): nsense_list3 = module.return_sense(nsense_list3, m, n)
						sense_list += nsense_list3

				for sense in sense_list:
					for j in xrange(0, len(sense)):
						if j == 0:
							fProb = float(fProb) * float(matrix_dic[j+1]['+']['+'+sense[j]])
							rProb = float(rProb) * float(matrix_dic[j+1]['-']['-'+sense[j]])
						else:
							fProb = float(fProb) * float(matrix_dic[j+1]['+'][sense[j-1]+sense[j]])
							rProb = float(rProb) * float(matrix_dic[j+1]['-'][sense[j-1]+sense[j]])

		for npPosReads in npReadPos_dic[npPos].keys():
			for npPosRead in npReadPos_dic[npPos][npPosReads]:
				if int(distance) <= 50: npPosRead.setreadRatio(forwardReads, reverseReads)
				else: npPosRead.setreadRatio(1, 1)
				if type == 'hmm': npPosRead.setreadProb(fProb, rProb)

	del npReadPos_dic, spReadPos_dic, npReadPos, spReadPos
	mediateMappingReads_dic = copy.deepcopy(npMappingReads_dic)
	for p in npMappingReads_dic.keys():
		readNum = len(npMappingReads_dic[p])
		if readNum > 1: #the mate is mapped
			matchpos_dic = dict()
			for m in range(len(npMappingReads_dic[p])):
				for n in range(len(npMappingReads_dic[p])):
					if (m < n and m != n):
						pairRead1 = npMappingReads_dic[p][m]; pairRead2 = npMappingReads_dic[p][n]
						if pairRead1.pos() == pairRead2.pnext() and pairRead1.pnext() == pairRead2.pos():
							matchpos_dic[m] = ''; matchpos_dic[n] = ''; reconfirm = 0
							
							if (float(pairRead1.readRatio()) == 1 and float(pairRead2.readRatio()) != 0) or (float(pairRead1.readRatio()) != 0 and float(pairRead2.readRatio()) == 1): reconfirm = 1	
							elif (float(pairRead1.readRatio()) == 0 and float(pairRead2.readRatio()) != 1) or (float(pairRead1.readRatio()) != 1 and float(pairRead2.readRatio()) == 0): reconfirm = -1
							else:
								if type == 'hmm':
									if (pairRead1.readProb() == 1 and pairRead2.readProb() == 1): pass
									elif (pairRead1.readProb() == 1111 and pairRead2.readProb() != -1111) or (pairRead1.readProb() != -1111 and pairRead2.readProb() == 1111): reconfirm = 1
									elif (pairRead1.readProb() == -1111 and pairRead2.readProb() != 1111) or (pairRead1.readProb() != 1111 and pairRead2.readProb() == -1111): reconfirm = -1
									elif (pairRead1.readProb() >= 1 and pairRead2.readProb() > 1) or (pairRead1.readProb() > 1 and pairRead2.readProb() >= 1): reconfirm = 1
									elif (pairRead1.readProb() <= 1 and pairRead2.readProb() < 1) or (pairRead1.readProb() < 1 and pairRead2.readProb() <= 1): reconfirm = -1
									else: pass

								if reconfirm == 0:
									if (float(pairRead1.readRatio()) >= 0.99 and float(pairRead2.readRatio()) >= 0.5) or (float(pairRead1.readRatio()) >= 0.5 and float(pairRead2.readRatio()) >= 0.99): reconfirm = 1
									elif (float(pairRead1.readRatio()) <= 0.01 and float(pairRead2.readRatio()) <= 0.5) or (float(pairRead1.readRatio()) <= 0.5 and float(pairRead2.readRatio()) <= 0.01): reconfirm = -1
									else:
										if (float(pairRead1.readRatio()) > 0.5 and float(pairRead2.readRatio()) >= 0.5) or (float(pairRead1.readRatio()) >= 0.5 and float(pairRead2.readRatio()) > 0.5): reconfirm = 1
										elif (float(pairRead1.readRatio()) < 0.5 and float(pairRead2.readRatio()) <= 0.5) or (float(pairRead1.readRatio()) <= 0.5 and float(pairRead2.readRatio()) < 0.5): reconfirm = -1
										else:
											if (float(pairRead1.readRatio()) == 0.5 and float(pairRead2.readRatio()) == 0.5): pass
											elif float(pairRead1.readRatio()) > 0.5:
												if 1-float(pairRead1.readRatio()) < float(pairRead2.readRatio()): reconfirm = 1
												else: reconfirm = -1
											elif float(pairRead1.readRatio()) < 0.5:
												if float(pairRead1.readRatio()) < 1-float(pairRead2.readRatio()): reconfirm = -1
												else: reconfirm = 1
											else: pass

							if reconfirm != 0:
								if reconfirm == 1:
									pairRead1.setreadRatio(1,0); pairRead2.setreadRatio(1,0)
								else: #reconfirm == -1:
									pairRead1.setreadRatio(0,1); pairRead2.setreadRatio(0,1)
								npairRead1 = module.flagCorrection(pairRead1); npairRead2 = module.flagCorrection(pairRead2)
								if not finalMappingReads_dic.has_key(p): finalMappingReads_dic[p] = []
								finalMappingReads_dic[p].extend([npairRead1, npairRead2])
							else: #reconfirm == 0:
								if not remainMappingReads_dic.has_key(p): remainMappingReads_dic[p] = []
								remainMappingReads_dic[p].extend([pairRead1, pairRead2])

			if len(matchpos_dic.keys()) > 0:
				matchposList = sorted(matchpos_dic.keys()); matchposList.reverse()
				if len(matchposList) == len(mediateMappingReads_dic[p]):
					del mediateMappingReads_dic[p]
				else:
					for matchpos in matchposList:
						del mediateMappingReads_dic[p][matchpos]

	del npMappingReads_dic
	for h in mediateMappingReads_dic.keys(): #the mate is unmapped
		for l in range(len(mediateMappingReads_dic[h])):
			singleRead = mediateMappingReads_dic[h][l]
			jpair = -1
			if singleRead.qname() in npExonJunctionReads_dic.keys():
				for s in range(len(npExonJunctionReads_dic[singleRead.qname()])):
					npjread = npExonJunctionReads_dic[singleRead.qname()][s]
					if singleRead.pos() == npjread.pnext() and singleRead.pnext() == npjread.pos():
						if 'XS:A:+' in npjread.tag(): singleRead.setreadRatio(1,0) #forward strand
						else: singleRead.setreadRatio(0,1) #reverse strand
						nsingleRead = module.flagCorrection(singleRead)
						if not npExonJunctionPairs_dic.has_key(singleRead.qname()):
							npExonJunctionPairs_dic[singleRead.qname()] = []
						npExonJunctionPairs_dic[singleRead.qname()].append(nsingleRead)
						jpair = 1; break

			if jpair < 0:
				jsingle = 0
				if float(singleRead.readRatio()) == 1: jsingle = 1
				elif float(singleRead.readRatio()) == 0: jsingle = -1
				else:
					if type == 'hmm':
						if singleRead.readProb() > 1: jsingle = 1
						elif singleRead.readProb() < 1: jsingle = -1
						else: pass
					if jsingle == 0:
						if float(singleRead.readRatio()) > 0.5: jsingle = 1
						elif float(singleRead.readRatio()) < 0.5: jsingle = -1
						else: pass

				if jsingle != 0:
					if jsingle == 1:
						singleRead.setreadRatio(1, 0)
					else: #jsingle == -1:
						singleRead.setreadRatio(0, 1)
					nsingleRead = module.flagCorrection(singleRead)
					if not finalMappingReads_dic.has_key(nsingleRead.qname()):
						finalMappingReads_dic[nsingleRead.qname()] = []
					finalMappingReads_dic[nsingleRead.qname()].append(nsingleRead)
				else: #jsingle == 0:
					if not remainMappingReads_dic.has_key(singleRead.qname()):
						remainMappingReads_dic[singleRead.qname()] = []
					remainMappingReads_dic[singleRead.qname()].append(singleRead)

	del mediateMappingReads_dic
	outputFile = open(outputdir + '/logs/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + '.logs', 'w')
	coutputSam = open(outputdir + '/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'c.sam', 'w')
	routputSam = open(outputdir + '/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'u.sam', 'w')
	module.writeSamHeader(coutputSam, assembly); module.writeSamHeader(routputSam, assembly)

	correctedReadsNum = 0; uncorrectedReadsNum = 0
	for a1 in finalMappingReads_dic.keys():
		for a2 in finalMappingReads_dic[a1]:
			coutputSam.write(module.writeSam(a2) + '\n')
			correctedReadsNum += 1
	for b1 in npExonJunctionReads_dic.keys():
		for b2 in npExonJunctionReads_dic[b1]:
			coutputSam.write(module.writeSam(b2) + '\n')
			correctedReadsNum += 1
	for c1 in npExonJunctionPairs_dic.keys():
		for c2 in npExonJunctionPairs_dic[c1]:
			coutputSam.write(module.writeSam(c2) + '\n')
			correctedReadsNum += 1
	for d1 in remainMappingReads_dic.keys():
		for d2 in remainMappingReads_dic[d1]:
			routputSam.write(module.writeSam(d2) + '\n')
			uncorrectedReadsNum += 1
	
	outputFile.write(chr + '\t' + str(correctedReadsNum) + '\t' + str(uncorrectedReadsNum) + '\n')
	outputFile.close(); coutputSam.close(); routputSam.close()

rpds(unstranded, stranded, chr)
sam.samToBam(outputdir + '/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'c.sam')
sam.bamSort(outputdir + '/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'c.bam')
sam.bamIndex(outputdir + '/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'c.bam')
