import requests

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

while True:

    response = requests.get ('http://localhost:7001/api/status/%s/' % jobid)
    status = response.json () ['status']

    print status
    if status not in ["PENDING", "STARTED", "RUNNING"]:
        break

if status != "SUCCESS":
    raise Exception (response.message)
