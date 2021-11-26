import re
from pprint import pprint

from hgvs.edit import Dup, Inv, NARefAlt
from hgvs.parser import Parser

RE_HGVS_G = r'[1-9]{1,2}:g\.(=|_|con|copy|del|dup|ins|inv|[0-9])'

# maps between chromosomes and refseq chromosome-level accessions
AC_MAP = {
    '1': 'NC_000001.10',
    '2': 'NC_000002.11',
    '3': 'NC_000003.11',
    '4': 'NC_000004.11',
    '5': 'NC_000005.9',
    '6': 'NC_000006.11',
    '7': 'NC_000007.13',
    '8': 'NC_000008.10',
    '9': 'NC_000009.11',
    '10': 'NC_000010.10',
    '11': 'NC_000011.9',
    '12': 'NC_000012.11',
    '13': 'NC_000013.10',
    '14': 'NC_000014.8',
    '15': 'NC_000015.9',
    '16': 'NC_000016.9',
    '17': 'NC_000017.10',
    '18': 'NC_000018.9',
    '19': 'NC_000019.9',
    '20': 'NC_000020.10',
    '21': 'NC_000021.8',
    '22': 'NC_000022.10',
    'X': 'NC_000023.10',
    '23': 'NC_000023.10',
    'Y': 'NC_000024.9',
}


class Hgvs:
    @staticmethod
    def __refseq_g_accession(hgvs_g_str):
        chromosome, rest = hgvs_g_str.split(":")
        return "%s:%s" % (AC_MAP[str(chromosome)], rest)

    @classmethod
    def parse(cls, hgvs_str):
        if re.match(RE_HGVS_G, hgvs_str):
            hgvsparser = Parser()
            ret = hgvsparser.parse_g_variant(
                cls.__refseq_g_accession(hgvs_str))

            print('--ret--')
            print(ret)

            print('--posedit--')
            pprint(ret.posedit)

            print('--pos--')
            pprint(ret.posedit.pos)

            print(f'--edit ({type(ret.posedit.edit)})--')
            pprint(ret.posedit.edit)

            if isinstance(ret.posedit.edit, NARefAlt):
                return {
                    'chromosome': hgvs_str.split(':')[0],
                    'start': ret.posedit.pos.start.base,
                    'end': ret.posedit.pos.end.base,
                    'ref': ret.posedit.edit.ref,
                    'alt': ret.posedit.edit.alt
                }
            elif isinstance(ret.posedit.edit, Dup):
                return {
                    'chromosome': hgvs_str.split(':')[0],
                    'start': ret.posedit.pos.start.base,
                    'end': ret.posedit.pos.end.base,
                    'ref': ret.posedit.edit.ref,
                    'alt': None
                }
            elif isinstance(ret.posedit.edit, Inv):
                return {
                    'chromosome': hgvs_str.split(':')[0],
                    'start': ret.posedit.pos.start.base,
                    'end': ret.posedit.pos.end.base,
                    'ref': None,
                    'alt': None
                }
