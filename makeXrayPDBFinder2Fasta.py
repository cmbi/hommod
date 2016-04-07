#!/usr/bin/python

import sys,os,urllib,re

scriptdir=os.path.abspath(os.path.dirname(sys.argv[0]))

sys.path.append(os.path.dirname(scriptdir)) # parent directory has modelutils.py

def writeFasta(seqs,path):

    f=open(path,'w')
    for key in seqs.keys():
        f.write('>%s\n%s\n'%(key,seqs[key]))
    f.close()

if len(sys.argv)!=2:
    print 'Usage: %s [output file]'%sys.argv[0]
    sys.exit(0)

outFile=sys.argv[1]

# Parse entire set of available sequences:
source=urllib.urlopen('ftp://ftp.cmbi.ru.nl/pub/molbio/data/pdbfinder2/PDBFIND2.TXT')
#source=open('/data/pdbfinder2/PDBFIND2.TXT','r')

# First gather all sequences from the pdbFinder:
pNuc=re.compile(r'[actgu\-]+')
pProt=re.compile(r'[ABCDEFGHIJKLMNOPQRSTUVWXYZ\-]+')

# Make a list of regular expressions to recognize unallowed HET-Groups
reUnallowedHetGroups=[]
for hetGroupName in open(os.path.join(scriptdir,'unallowed-hetgroups'),'r'):

    hetGroupName=hetGroupName[:-1] # Remove the newline character
    if len(hetGroupName.strip())>0:
        pattern=re.compile(hetGroupName)
        reUnallowedHetGroups.append(pattern)

errorTemplates=open(os.path.join(scriptdir,'blacklisted_templates'),'r').read().lower().split()

seqs={}
hetatoms={}

chainID=None
pdbID=None
section=None
expMethod=None
resolution=None

line = source.readline()
while len(line) > 0:

    s=line.split(':')

    if not line[0].isspace():
        section=s[0].strip() # Chain, HET-Groups, Exp-Method, etc.

    if s[0].strip()=='ID':
        pdbID=s[1].strip().lower()

        if pdbID in errorTemplates:

            # Go to next entry:
            while line.strip()!='//':
                line=source.readline()

    elif s[0].strip()=='Chain':

        chainID=s[1].strip()

    elif s[0].strip()=='Exp-Method':

        expMethod=s[1].strip()

        if expMethod!='X':

            # Go to next entry:
            while line.strip () != '//' and len (line) > 0:
                line=source.readline()

            expMethod=None

    elif s[0].strip()=='Resolution' and section=='Exp-Method':

        resolution = float(s[1])

    elif s[0].strip()=='Name' and section=='HET-Groups':

        hetgroupname=s[1].strip()

        for pattern in reUnallowedHetGroups:

            m=pattern.search( hetgroupname )
            if m:
                # Go to next entry:
                while line.strip() != '//' and len (line) > 0:
                    line=source.readline()

                break

    elif s[0].strip()=='Sequence':

        seq=s[1].strip()

        if len(seq)>0 and expMethod=='X' and pProt.match(seq) and not pNuc.match(seq) :

            # Generate key, unique for pdb file and chain id:
            key='pdb|%s|%s'%(pdbID.upper(),chainID)
#            key='%s:%s|PDBID|CHAIN|SEQUENCE'%(pdbID.upper(),chainID)

            # Take the longest sequence in the chain:
            if not seqs.has_key(key) or len(seqs[key])<len(seq):

                seqs[key] = seq.replace('-','X')

    line=source.readline()

source.close()

writeFasta (seqs, outFile)
