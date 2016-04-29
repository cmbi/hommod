from hommod_rest.factory import create_app, create_celery_app
from time import sleep

import json

class TestTasks(object):

    @classmethod
    def setupClass(self):

        import hommod_rest.default_settings as config

        self.app = create_app({'TESTING': True, \
                               'CELERY_ALWAYS_EAGER': True,
                               'YASARADIR': config.YASARADIR})
        self.celery = create_celery_app(self.app)

    def test_job(self):

        response = self.app.test_client().post('/api/submit/',data={'position': 100,
                                                              'species_id': "HUMAN",
                                                              'sequence': \
    "MGAARGSPARPRRLPLLSVLLLPLLGGTQTAIVFIKQPSSQDALQGRRALLRCEVEAPGPVHVYWLLDGAPVQD" + \
    "TERRFAQGSSLSFAAVDRLQDSGTFQCVARDDVTGEEARSANASFNIKWIEAGPVVLKHPASEAEIQPQTQVTL" + \
    "RCHIDGHPRPTYQWFRDGTPLSDGQSNHTVSSKERNLTLRPAGPEHSGLYSCCAHSAFGQACSSQNFTLSIADE" + \
    "SFARVVLAPQDVVVARYEEAMFHCQFSAQPPPSLQWLFEDETPITNRSRPPHLRRATVFANGSLLLTQVRPRNA" + \
    "GIYRCIGQGQRGPPIILEATLHLAEIEDMPLFEPRVFTAGSEERVTCLPPKGLPEPSVWWEHAGVRLPTHGRVY" + \
    "QKGHELVLANIAESDAGVYTCHAANLAGQRRQDVNITVATVPSWLKKPQDSQLEEGKPGYLDCLTQATPKPTVV" + \
    "WYRNQMLISEDSRFEVFKNGTLRINSVEVYDGTWYRCMSSTPAGSIEAQARVQVLEKLKFTPPPQPQQCMEFDK" + \
    "EATVPCSATGREKPTIKWERADGSSLPEWVTDNAGTLHFARVTRDDAGNYTCIASNGPQGQIRAHVQLTVAVFI" + \
    "TFKVEPERTTVYQGHTALLQCEAQGDPKPLIQWKDKPVPEESEGPGSPPPYKMIQTIGLSVGAAVAYIIAVLGL" + \
    "MFYCKKRCKAKRLQKQPEGEEPEMECLNGGPLQNGQPSAEIQEEVALTSLGSGPAATNKRHSTSDKMHFPRSSL" + \
    "QPITTLGKSEFGEVFLAKAQGLEEGVAETLVLVKSLQSKDEQQQLDFRRELEMFGKLNHANVVRLLGLCREAEP" + \
    "HYMVLEYVDLGDLKQFLRISKSKDEKLKSQPLSTKQKVALCTQVALGMEHLSNNRFVHKDLAARNCLVSAQRQV" + \
    "KVSALGLSKDVYNSEYYHFRQAWVPLRWMSPEAILEGDFSTKSDVWAFGVLMWEVFTHGEMPHGGQADDEVLAD" + \
    "LQAGKARLPQPEGCPSKLYRLMQRCWALSPKDRPSFSEIASALGDSTVDSKP" })

        jobid = json.loads (response.get_data ()) ['jobid']

        prevstatus = ''
        while True:

            response = self.app.test_client().get ('api/status/%s/'% jobid)
            status = json.loads (response.get_data ()) ['status']

            if status != prevstatus:
                print status

            prevstatus = status

            if status not in ["PENDING", "STARTED", "RUNNING"]:
                break

            sleep (10)

        if status != "SUCCESS":

            js = json.loads (response.get_data ())
            if "message" in js:
                raise Exception (js ['message'])
            else:
                raise Exception ("status is " + status)

        response = self.app.test_client ().get ('api/get_model_file/%s.pdb' % jobid)
        response = self.app.test_client ().get ('api/get_metadata/%s/' % jobid)

