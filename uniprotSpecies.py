#!/usr/bin/python

import sys,os,subprocess

from gzip import GzipFile
from urllib2 import urlopen
from io import BufferedReader

fastas = ['trembl','sprot']

if len(sys.argv) != 2:
	print 'Usage: %s <directory>'%sys.argv[0]
	sys.exit(0)

outDir = sys.argv[1]
if not os.path.isdir(outDir):
	os.mkdir(outDir)

pid = os.getpid()

tmpfiles={}
def tmpStore(ID,seq,species):

	path = '/tmp/uniprot-%i-%s.fa'%(pid,species)
	open( path , 'a').write( '>%s\n%s\n'%(ID,seq) )
	tmpfiles[species] = path

for db in fastas:
	url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_%s.fasta.gz' % db

	currentID = None
	currentSeq= None
	currentSpecies = None

	proc = subprocess.Popen('curl -s %s | zcat' % url,shell=True,stdout=subprocess.PIPE)

	for line in proc.stdout:

		if line.startswith(">"):

			if currentSeq:
				tmpStore(currentID,currentSeq,currentSpecies)

			currentID = line[1:].split()[0]
			currentSeq = ''
			currentSpecies = currentID.split('_')[1]
		elif len(line.strip())>0:
			currentSeq += line.strip()

	if currentSeq:
		tmpStore(currentID,currentSeq,currentSpecies)

for species in tmpfiles.keys():
	outFasta = os.path.join(outDir,'uniprot-%s.fa'%species)
	os.system('mv %s %s'%(tmpfiles[species],outFasta))
