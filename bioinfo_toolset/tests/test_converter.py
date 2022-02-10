from unittest import TestCase

from bioinfo_toolset.modules.converter import three_to_one


class TestThreeToOne(TestCase):
    def test_three_to_one(self):
        cases = {
            'p.Val600Glu': 'p.V600E',
            'p.Asp698GlufsTer22': 'p.D698EfsTer22'
        }

        for three, expected in cases.items():
            self.assertEqual(expected, three_to_one(three))
