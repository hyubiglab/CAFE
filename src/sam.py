#!/share/apps/python/bin/python
import sys, os
from subprocess import Popen, PIPE, STDOUT

def samToBam(samfile):
	bamfile = '.'.join(samfile.split('.')[:-1]) + '.bam'
	commandBase = 'samtools view -Sb ' + samfile + ' > ' + bamfile
	print commandBase
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()
	os.remove(samfile)

def bamRetrieve(bamfile, locus):
	commandBase = 'samtools view ' + bamfile + ' ' + locus
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	output=p.stdout.read()
	return output

def bamSort(bamfile):
	sortedfile = '.'.join(bamfile.split('.')[:-1]) + 'sorted'
	commandBase = 'samtools sort -@ 4 ' + bamfile + ' ' + sortedfile ##
	print commandBase
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()
	print bamfile
	os.remove(bamfile)
	print sortedfile + '.bam', bamfile
	os.renames(sortedfile + '.bam', bamfile)

def bamIndex(bamfile):
	commandBase = 'samtools index ' + bamfile
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def bamCat(bamfiles, bamfile):
	commandBase = 'samtools cat -o ' + bamfile + ' ' + ' '.join(bamfiles)
	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()
	bamSort(bamfile)
	bamIndex(bamfile)
