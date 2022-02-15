from unittest import TestCase
from bioinfo_toolset.modules.vep import vep
from bioinfo_toolset.modules.formatter import transcript_name


class TestFormatter(TestCase):
    def test_transcript_name(self):
        """Test conversion of transcriptnames"""
        transcript_map_chgr37 = {
            '5:149439278-149439301/GC': 'D698EfsTer22',
            '17:7577610-7577610/C': 'splice site 673-2A>G'
        }

        for region, expected in transcript_map_chgr37.items():
            vep_output = vep(region, input_type='region', GRCh37=True)
            for transcript in vep_output[0]['transcript_consequences']:
                # We only check canonical transcripts
                if transcript.get('canonical'):
                    self.assertEqual(expected, transcript_name(transcript))
