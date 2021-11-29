from unittest import TestCase

from bioinfo_toolset.modules.converter import three_to_one

class TestThreeToOne(TestCase):
    def test_three_to_one(self):
        self.assertEqual('p.V600E', three_to_one('p.Val600Glu'))