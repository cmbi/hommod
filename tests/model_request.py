import requests
from time import sleep

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
    raise Exception (response.message)

response = requests.get ('http://localhost:7001/api/get_model_file/%s.pdb' % jobid)
print response.text
