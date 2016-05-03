#! /usr/bin/python

# This script can be used to submit a model request to a
# running service. It's made for testing purposes.

import logging
import sys
import os
import traceback
import requests

_log = logging.getLogger(__name__)

if len(sys.argv) != 4:
    print 'usage: %s [sequence] [species id] [position]' % sys.argv[0]
    sys.exit(1)

species = sys.argv[2].upper()
seq = sys.argv[1].upper()
pos = int (sys.argv[3])

# Send http post request to server:

url = 'http://localhost:7001/api'

payload = {'sequence': seq, 'species_id': species, 'position': pos}
response = requests.post (url + '/submit/', data = payload)
