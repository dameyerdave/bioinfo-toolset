from pprint import pprint
from unittest import TestCase
from bioinfo_toolset.modules.hgvs import Hgvs


class TestHgvs(TestCase):
    def setUp(self):
        self.hgvs = Hgvs()

    def test_from_transcript_change(self):
        hgvs_map_grch37 = {
            '5:149439278:2094_2117>GC': '5:149439278-149439301/GC',
            'X:53239699:1642_1643GT>AA': 'X:53239699-53239700/TT',
            '8:37555623:1204_1205insCGACGCACCTCGGCGGCTCCAGCTGCTCCACCTGCAGCGCGC': '8:37555624-37555623/CGACGCACCTCGGCGGCTCCAGCTGCTCCACCTGCAGCGCGC',
            'X:129190010:5035_5036insT': 'X:129190011-129190010/T',
            '17:7577610:673-2A>G': '17:7577610-7577610/C',
        }

        for test_transcript, expected_region in hgvs_map_grch37.items():
            result = self.hgvs.from_transcript_change(
                *test_transcript.split(':'), True)

            self.assertEqual(
                expected_region, result['region'], 'region is not correct')
