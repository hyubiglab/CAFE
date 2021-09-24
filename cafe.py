#!/share/apps/python/bin/python
import sys, os
from subprocess import Popen, PIPE, STDOUT

import config as conf

if len(sys.argv) < 7:
	print '\nCAFE.1.0.1: Co-Assembly of stranded and unstranded RNA-seq data Followed by End-correction\n'
	print 'usage: python cafe.py [options] [parameters]...'
	print 'ex) python cafe.py -g human -v hg19 -c all -t npsp -u unstranded.bam -s stranded.bam -o outputdir\n'
	print 'options:'
	print '%-24s %-50s'%(' -g <genome>','genome name (human or mouse)')
	print '%-24s %-50s'%(' -v <assembly>','assembly version (hg19 or mm9)')
	print '%-24s %-50s'%(' -c <chr>','chromosome to assembly (all or specific chromosome)' )
	print '%-24s %-50s'%(' -t <type>','type of input (npsp, np, or sp)')
	print '%-24s %-50s'%(' -u <unstranded.bam>','unstranded RNA-seq bamfile')
	print '%-24s %-50s'%(' -s <stranded.bam>','stranded RNA-seq bamfile')
	print '%-24s %-50s'%(' -o <outputdir>','output directory')
	print '%-24s %-50s'%(' -h, --help','show this screen\n')
	sys.exit(0)

for i in xrange(len(sys.argv)):
	if sys.argv[i]=='-g': genome = sys.argv[i + 1]
	elif sys.argv[i]=='-v': assembly = sys.argv[i+ 1]
	elif sys.argv[i]=='-c': chrs = sys.argv[i+ 1]
	elif sys.argv[i]=='-t': type = sys.argv[i + 1]
	elif sys.argv[i]=='-u': unstranded = sys.argv[i + 1]
	elif sys.argv[i]=='-s': stranded = sys.argv[i + 1]
	elif sys.argv[i]=='-o': outputdir = sys.argv[i + 1]
	elif sys.argv[i]=='-h' or sys.argv[i]=='--help':
		command = 'python cafe.py'
                p=Popen(command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
                p.wait()
                sys.exit(0)

if chrs == 'all':
	if assembly == 'hg19': chrs = conf.hm_chrs
	elif assembly == 'mm9': chrs = conf.ms_chrs
else: chrs = [chrs]

def makedir(outputdir, dir):
	noutputdir = outputdir + '/' + dir + '/'; nlogdir = noutputdir + 'logs/'
	if not os.path.exists(nlogdir): os.makedirs(nlogdir)
	return noutputdir

def execute(command):
	print command
	p=Popen(command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def cufflinks(assembly, bamfile, outputdir, libType = 'fr-secondstrand'):
	if assembly == 'hg19': intronmin = conf.hm_intronmin; intronmax = conf.hm_intronmax
	elif assembly == 'mm9': intronmin = conf.ms_intronmin; intronmax = conf.ms_intronmax
	command = 'cufflinks -p 4 --min-intron-length ' + str(intronmin) + ' --max-intron-length ' + str(intronmax) + ' --library-type ' + libType + ' --no-update-check -o ' + outputdir + ' ' + bamfile
	print command
	p=Popen(command, shell=True)
	p.wait()

workingdir = os.getcwd()

# rpds
noutputdir = makedir(outputdir, 'rpds')
for chr in chrs: execute('python ' + workingdir + '/src/rpds.py ' + assembly + ' ' + unstranded + ' ' + stranded + ' ' + chr + ' ' + noutputdir)

# assembly
noutputdir = makedir(outputdir, 'gtf/rpds/') #rpds
for chr in chrs: cufflinks(assembly, outputdir + '/rpds/' + unstranded.split('/')[-1][:-4].split('_')[0] + '_' + chr + 'c.bam', noutputdir)

noutputdir = makedir(outputdir, 'gtf/stranded/') #stranded
cufflinks(assembly, stranded, noutputdir)

# combine
noutputdir = makedir(outputdir, 'gtf/combine/')
for chr in chrs: execute('python ' + workingdir + '/src/combine.py ' + assembly + ' ' + outputdir + '/gtf/rpds/transcripts.gtf,' + outputdir + '/gtf/stranded/transcripts.gtf ' + chr + ' ' + noutputdir)

# exj
noutputdir = makedir(outputdir, 'gtf/exj/')
for chr in chrs: execute('python ' + workingdir + '/src/exj.py ' + assembly + ' ' + outputdir + '/gtf/combine/transcripts_' + chr + '.gtf ' + outputdir + '/rpds/' + unstranded.split('/')[-1][:-4] + 'c.bam ' + stranded + ' ' + chr + ' ' + noutputdir)

# tss
noutputdir = makedir(outputdir, 'gtf/tss/')
for chr in chrs: execute('python ' + workingdir + '/src/tss.py ' + assembly + ' ' + outputdir + '/gtf/exj/transcripts_' + chr + '.exj.gtf ' + chr + ' ' + noutputdir)

# cps
noutputdir = makedir(outputdir, 'gtf/cps/')
for chr in chrs: execute('python ' + workingdir + '/src/cps.py ' + assembly + ' ' + outputdir + '/gtf/tss/transcripts_' + chr + '.tss.gtf ' + chr + ' ' + noutputdir)

# filter
noutputdir = makedir(outputdir, 'gtf/filter/')
for chr in chrs: execute('python ' + workingdir + '/src/filter.py ' + assembly + ' ' + outputdir + '/gtf/cps/transcripts_' + chr + '.tss.cps.gtf ' + chr + ' ' + noutputdir)
