from nose.tools import eq_

from hommod.services.interpro import interpro



def test_interpro_parser():

    xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<protein-matches xmlns="http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5">
    <protein>
        <sequence md5="?">FAKE</sequence>
        <matches>
        </matches>
    </protein>
</protein-matches>"""

    domain_ranges = interpro._parse_interpro_ranges(xml)
    eq_(type(domain_ranges), list)
