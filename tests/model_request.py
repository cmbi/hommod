import requests

payload = {
    'position': 12, 'species_id': 'RAT',
    'sequence': """
MLATRALSLIGKRAISTSVCLRAHGSVVKSEDYALPSYVDRRDYPLPDVAHVKLLSASQK
ALKEKEKADWSSLSRDEKVQLYRIQFNESFAEMNKGTNEWKTVVGLAMFFIGFTALVLIW
EKSYVYGPIPHTFDRDWVAMQTKRMLDMKVNPIQGFSAKWDYNKNEWKK
""".replace('\n', '')
}

response = requests.post('http://localhost:6000/api/submit/', data=payload)

print response.text
