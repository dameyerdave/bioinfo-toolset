from unittest import TestCase
from bioinfo_toolset.modules.helper import inverse


class Test_Helper(TestCase):
    def test_inverse(self):
        """Test inverse a string"""
        self.assertEqual('GTCA', inverse('ACTG'))
