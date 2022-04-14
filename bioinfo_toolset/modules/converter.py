from Bio.PDB.Polypeptide import three_to_one as biopython_three_to_one, one_to_three as biopython_one_to_three, aa3
import re


def three_to_one(three_letter_string: str):
    converted = three_letter_string
    for aa_name in aa3:
        converted = re.sub(aa_name, biopython_three_to_one(
            aa_name), converted, flags=re.I)
    return converted


def one_to_three(one_letter_string: str):
    match = re.search(
        r'(?P<prefix>p.)?(?P<from>[ACDEFGHIKLMNPQRSTVWY])(?P<codon>[0-9]+)(?P<to>[ACDEFGHIKLMNPQRSTVWY])(?P<postfix>.*)', one_letter_string)
    if match:
        _from = biopython_one_to_three(match.group('from')).capitalize()
        _to = biopython_one_to_three(match.group('to')).capitalize()
        return f"{match.group('prefix') if match.group('prefix') else ''}{_from}{match.group('codon')}{_to}{match.group('postfix') if match.group('postfix') else ''}"
    return one_letter_string
