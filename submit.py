#! /usr/bin/python

# This script can be used to submit a model request to a
# running service. It's made for testing purposes.

import logging
import sys
import os
import traceback
import requests
from hommod_rest.services.modelutils import parseFirstFastaSequence

_log = logging.getLogger(__name__)

if len(sys.argv) != 4 and len(sys.argv) != 5:
    print 'usage: %s [sequence/fasta path] [species id] [position] [template id]' % sys.argv[0]
    sys.exit(1)

species = sys.argv[2].upper()

fasta = sys.argv[1]
if os.path.isfile (fasta):
    _id, seq = parseFirstFastaSequence (fasta)
else:
    seq = fasta.upper()

pos = int (sys.argv[3])

# Send http post request to server:

url = 'http://localhost:7001/api'

payload = {'sequence': seq, 'species_id': species, 'position': pos}
if len(sys.argv) != 5:
    payload['template_id'] = sys.argv[4]

response = requests.post (url + '/submit/', data = payload)

print response.json () ['jobid']
