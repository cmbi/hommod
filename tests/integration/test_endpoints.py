from time import sleep

import json

import logging
_log = logging.getLogger (__name__)
_log.setLevel (logging.DEBUG)

import requests

url = "http://localhost:7001/api"

def _test_func (sequence, species_id, position):

    _log.debug ("submitting..")

    response = requests.post ('%s/submit/' % url,
                              data = {'position'  : position,
                                      'species_id': species_id,
                                      'sequence'  : sequence})

    jobid = response.json () ['jobid']

    _log.debug ("recieved job id: " + jobid)

    prev_status = ''
    while True:

        response = requests.get ('%s/status/%s/'% (url, jobid))
        status = response.json () ['status']

        if status != prev_status:
            _log.debug ("got status %s for job %s " %(status, jobid))
        prev_status = status

        if status in ["PENDING", "STARTED", "RUNNING"]:

            sleep (1)
        else:
            break


    if status != "SUCCESS":

        js = response.json ()
        if "message" in js:
            raise Exception (js ['message'])
        else:
            raise Exception ("job %s status is %s" % (jobid, status))

    response = requests.get ('%s/get_model_file/%s.pdb' % (url, jobid))
    response = requests.get ('%s/get_metadata/%s/' % (url, jobid))


# Best use very easy requests, that don't take much time!

def test_endpoints_1crn ():

    _test_func ("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", "CRAAB", 25)
