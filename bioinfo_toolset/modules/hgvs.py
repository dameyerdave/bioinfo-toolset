import re
import traceback

from hgvs.edit import Dup, Inv, NARefAlt
from hgvs.parser import Parser

from bioinfo_toolset.modules.vrs import VRS
from bioinfo_toolset.modules.lookup import complement_allele_lookup

from friendlylog import colored_logger as log

from bioinfo_toolset.modules.helper import inverse

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

RE_HGVS_G = r'([1-9]{1}(?:[0-9]{1})?|[XY]{1}):g\.(=|_|con|copy|del|dup|ins|inv|[0-9])'
RE_HGVS_C = r'(?:c\.)?(?P<position>[^ACTG]+)(?P<from_allele>[ACTG]+)>(?P<to_allele>[ACTG]+)'
RE_TRANS_C = [
    r'(?:c\.)?(?P<position>[^ACTG]+)(?P<from_allele>[ACTG]+)>(?P<to_allele>[ACTG]+)',
    r'(?:c\.)?(?P<position>[0-9]+(?:_[0-9]+)?)(?P<type>(ins|delins|del|dub|inv|>))(?P<to_allele>[ACTG]+)',
    r'(?:c\.)?(?P<position>[0-9+-]+(?:_[0-9+-]+)?)(?P<type>(del|dup|inv))(?P<num_alleles>[0-9]*)'
]
RE_POS_DELIM = r'_'


class Hgvs:
    hgvsparser = Parser()
    vrs = VRS()

    @staticmethod
    def __refseq_g_accession(hgvs_g_str):
        chromosome, rest = hgvs_g_str.split(":")
        return "%s:%s" % (AC_MAP[str(chromosome)], rest)

    @classmethod
    def from_transcript_change(cls, chromosome: str, position: int, transcript_change: str, GRCh37: bool = False):
        for rex in RE_TRANS_C:
            if re.match(rex, transcript_change):
                try:
                    position = int(position)
                    transcript_change_info = re.search(rex, transcript_change)
                    if '_' in transcript_change_info.group('position'):
                        start, end = list(map(lambda p: int(eval(p)), re.split(RE_POS_DELIM, transcript_change_info.group(
                            'position'))))
                        position_part = f"{position}_{position + end - start}"
                    else:
                        position_part = position
                    if 'from_allele' in transcript_change_info.re.groupindex:
                        # we found a change nomenclature
                        reference_allele = cls.vrs.allele_at_position(
                            'GRCh37' if GRCh37 else 'GRCh38', chromosome, position - 1, position - 1 + len(transcript_change_info.group('from_allele')))
                        complement_allele = complement_allele_lookup(
                            reference_allele)
                        ref = reference_allele
                        if transcript_change_info.group('from_allele') == reference_allele:
                            alt = transcript_change_info.group(
                                'to_allele')
                        elif inverse(transcript_change_info.group('from_allele')) == reference_allele:
                            alt = inverse(transcript_change_info.group(
                                'to_allele'))
                        # Sanity check it should then be the to_allele
                        elif transcript_change_info.group('from_allele') == complement_allele:
                            # we need to invert the values (because its on the backwards strand)
                            alt = complement_allele_lookup(transcript_change_info.group(
                                'to_allele'))
                        elif inverse(transcript_change_info.group('from_allele')) == complement_allele:
                            alt = inverse(complement_allele_lookup(transcript_change_info.group(
                                'to_allele')))
                        else:
                            log.warning(
                                f"Something must be wrong the from allele ({transcript_change_info.group('from_allele')}) does neather correspond to the reference allele ({reference_allele}) not to the complement allele ({complement_allele}): {chromosome}:{position} {transcript_change}")
                            ref = transcript_change_info.group('from_allele')
                            alt = transcript_change_info.group('to_allele')

                        return cls.parse(f"{chromosome}:g.{position_part}{ref}>{alt}")
                    elif 'type' in transcript_change_info.re.groupindex:
                        # we found a insertion, deletion...
                        if transcript_change_info.group('type') in ['ins', 'delins']:
                            # VEP requirement: An insertion (of any size) is indicated by start coordinate = end coordinate + 1
                            position_part = f"{position + 1}_{position}"
                            return cls.parse(f"{chromosome}:g.{position_part}{transcript_change_info.group('type')}{transcript_change_info.group('to_allele')}")
                        elif transcript_change_info.group('type') in ['del', 'dup', 'inv']:
                            if 'to_allele' in transcript_change_info.re.groupindex:
                                reference_allele = cls.vrs.allele_at_position(
                                    'GRCh37' if GRCh37 else 'GRCh38', chromosome, position, position + len(transcript_change_info.group('to_allele')))
                                complement_allele = complement_allele_lookup(inverse(cls.vrs.allele_at_position(
                                    'GRCh37' if GRCh37 else 'GRCh38', chromosome, position, position + len(transcript_change_info.group('to_allele')))))
                                if transcript_change_info.group('to_allele') == reference_allele:
                                    position_part = f"{position + 1}_{position + len(transcript_change_info.group('to_allele'))}"
                                elif transcript_change_info.group('to_allele') == complement_allele:
                                    position_part = f"{position + 1}_{position + len(transcript_change_info.group('to_allele'))}"
                            return cls.parse(f"{chromosome}:g.{position_part}{transcript_change_info.group('type')}")
                        elif transcript_change_info.group('type') in ['>']:
                            return cls.parse(f"{chromosome}:g.{position_part}delins{transcript_change_info.group('to_allele')}")
                except Exception as ex:
                    log.error(ex)
                    traceback.print_exc()
        log.error(
            f"Cannot get variant from {chromosome}:{position}:{transcript_change}: No regex matched.")
        return None

    @ classmethod
    def parse_c(cls, hgvs_str):
        """Parses a hgvs c string. A hgvs string at the transcript level."""
        # if re.match(RE_HGVS_C, hgvs_str):
        #     cls.hgvsparser.parse_hgvs_variant()
        raise NotImplementedError()

    @ classmethod
    def parse_g(cls, hgvs_str):
        """Same as parse()"""
        return cls.parse(hgvs_str)

    @ classmethod
    def parse(cls, hgvs_str):
        """Parses a hgvs g string. A hgvs string at the gene level."""
        if re.match(RE_HGVS_G, hgvs_str):
            try:
                ret = cls.hgvsparser.parse(
                    cls.__refseq_g_accession(hgvs_str))

                # if isinstance(ret.posedit.edit, NARefAlt):
                chromosome = hgvs_str.split(':')[0]
                position = re.sub(r'_', '-', str(ret.posedit.pos))
                allele = re.sub(r'>', '/', str(ret.posedit.edit).upper())
                return {
                    'chromosome': chromosome,
                    'start': ret.posedit.pos.start.base,
                    'end': ret.posedit.pos.end.base,
                    'ref': ret.posedit.edit.ref,
                    'alt': ret.posedit.edit.alt,
                    'region': f"{chromosome}:{ret.posedit.pos.start.base}-{ret.posedit.pos.end.base}/{ret.posedit.edit.alt if ret.posedit.edit.alt else allele}",
                    'identifier': f"{chromosome}_{ret.posedit.pos.start.base}_{ret.posedit.pos.end.base}{'_' + ret.posedit.edit.ref if ret.posedit.edit.ref else ''}_{'-/' + ret.posedit.edit.alt if ret.posedit.pos.start.base > ret.posedit.pos.end.base else ret.posedit.edit.alt if ret.posedit.edit.alt else allele}"
                }
                # elif isinstance(ret.posedit.edit, Dup):
                #     return {
                #         'chromosome': hgvs_str.split(':')[0],
                #         'start': ret.posedit.pos.start.base,
                #         'end': ret.posedit.pos.end.base,
                #         'ref': ret.posedit.edit.ref,
                #         'alt': None,
                #         'region': 'TODO',
                #         # 'vcf': 'DUP'
                #     }
                # elif isinstance(ret.posedit.edit, Inv):
                #     return {
                #         'chromosome': hgvs_str.split(':')[0],
                #         'start': ret.posedit.pos.start.base,
                #         'end': ret.posedit.pos.end.base,
                #         'ref': None,
                #         'alt': None,
                #         'region': 'TODO',
                #         # 'vcf': 'INV'
                #     }
            except Exception as ex:
                log.warning(f"HGVSParser error: {ex}")
                chr, rest = hgvs_str.split(':')
                rest = re.sub(r'g\.', '', rest)
                info = re.search(
                    r'(?P<position>[0-9_]+)(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)', rest)
                if '_' in info.group('position'):
                    start, end = re.split(RE_POS_DELIM, info.group('position'))
                else:
                    start = end = info.group('position')
                return {
                    'chromosome': chr,
                    'start': start,
                    'end': end,
                    'ref': info.group('ref'),
                    'alt': info.group('alt'),
                    'region': f"{chr}:{start}-{end}/{info.group('alt')}",
                    'identifier': f"{chr}_{start}_{end}_{info.group('ref')}_{info.group('alt')}"
                    # 'vcf': f"{chr} {start} {end} {info.group('ref')}/{info.group('alt')}"
                }
        # In cases we do not find a match we return None
        log.warning(f"No matching variant found for {hgvs_str}")
        return None
