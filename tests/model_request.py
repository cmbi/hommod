import requests
from time import sleep
import sys

payload = {
    'position': 47, 'species_id': 'HUMAN',
    'sequence': """
VHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIGNAKVK
AHGKKVITAFNDGLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNIVIVLGHHLGKDF
TPAAQAAFQKVVAGVATALAHKYH
""".replace('\n', '')
}

response = requests.post('http://localhost:7001/api/submit/', data=payload)

jobid = response.json () ['jobid']

prevstatus = ''
while True:

    response = requests.get ('http://localhost:7001/api/status/%s/' % jobid)
    status = response.json () ['status']

    if status != prevstatus:
        print status

    prevstatus = status

    if status not in ["PENDING", "STARTED", "RUNNING"]:
        break

    sleep (10)

if status != "SUCCESS":
    if "message" in response.json ():
        raise Exception (response.json () ['message'])
    else:
        sys.exit (1)

response = requests.get ('http://localhost:7001/api/get_model_file/%s.pdb' % jobid)
print response.text
