#!/usr/bin/env python

from bioinfo_toolset.modules.vep_offline import OfflineVep
import traceback
from pprint import pprint
from sys import stdout

import click
from click import Choice
from friendlylog import colored_logger as log
from termcolor import colored


@click.group()
def cli():
    pass


@click.command()
@click.option('--from', '-f', '_from', default='hg19', type=str, help='Convert from this genome assembly')
@click.option('--to', '-t', default='hg38', type=str, help='Convert to this genome assembly')
@click.argument('chromosome', type=int)
@click.argument('position', type=int)
def liftover(_from, to, chromosome, position):
    from bioinfo_toolset.modules.liftover import liftover
    try:
        t_chr, t_pos = liftover(_from, to, chromosome, position)
        source_str = "%2s:%-10s" % (chromosome, position)
        target_str = "%2s:%-10s" % (t_chr, t_pos)
        print("%-16s -> %-16s" % (_from, to))
        print("%-16s -> %-16s" % (source_str, target_str))
    except ValueError as ve:
        log.error(ve)
        traceback.print_exc(file=stdout)


cli.add_command(liftover)


@click.command()
@click.option('--parse', '-p', '_parse', is_flag=True, help='Parse the hgvs string')
@click.argument('hgvs_string', type=str)
def hgvs(hgvs_string, _parse):
    from bioinfo_toolset.modules.hgvs import Hgvs
    if _parse:
        pprint(Hgvs.parse(hgvs_string))


cli.add_command(hgvs)


@click.command(context_settings={"ignore_unknown_options": True}, help="The query string can be a hgvs, id or region (e.g. hgvs: 9:g.22125504G>C, region: 3:178928079/C)")
@click.option('--GRCh37', '--grch37', '--old', '-o', 'GRCh37', is_flag=True, help="Use the GRCh rest api for the query")
@click.option('--liftover', '-lift', '-l', '_liftover', is_flag=True, help="Liftover the variant GRCh38/GRCh37, depending on the used assembly")
@click.option('--species', '-s', default='human', type=str, help='The species to query VEP')
@click.option('--input-type', '-t', 'input_type', default='hgvs', type=Choice(['hgvs', 'id', 'region', 'vcf', 'vcf_file']), help='The input type you want to query')
@click.option('--enrich-transcripts', '-e', 'enrich_transcripts', is_flag=True, help='Enrich the transcripts (takes longer)')
@click.option('--all-transcripts', '-a', 'all_transcripts', is_flag=True, help='Only show canonical transcript(s)')
@click.option('--refseq-mode', '-r', 'refsec_mode', is_flag=True, help='Use RefSeq transcript set to report consequences')
@click.argument('input', type=str)
def vep(species, input_type, input, GRCh37, _liftover, enrich_transcripts, all_transcripts, refsec_mode):
    from bioinfo_toolset.modules.vep import vep
    from bioinfo_toolset.modules.liftover import liftover

    def output(text, indent=0):
        print(f'%-{indent * 4}s{text}' % (' '))

    def output_kv(id, value, indent=0, highlight=False):
        attrs = []
        if highlight:
            attrs = ['bold']
        if isinstance(value, list):
            value = ', '.join(value)
        print(f'%-{indent * 4}s%-{34 if highlight else 30}s : %s' %
              (' ', colored(id, 'cyan', attrs=attrs), colored(value, 'white', attrs=attrs)))

    def output_item(_dict, id, indent=0, highlight=False):
        if id in _dict:
            output_kv(id, _dict[id], indent, highlight)

    def output_splitter(indent=0):
        output('---', indent)

    def output_variant(variant, indent=0):
        for id in ['assembly_name',
                   'id',
                   'seq_region_name',
                   'start',
                   'end',
                   'allele_string',
                   'strand',
                   'somatic',
                   'most_severe_consequence',
                   'clin_sig',
                   'frequencies']:
            output_item(variant, id, indent)

    def output_transcript(transcript, indent=0):
        highlight = False
        if 'canonical' in transcript and transcript['canonical'] == 1:
            highlight = True

        for id in [
            'transcript_id',
            'gene_symbol',
            'gene_id',
            'variant_allele',
            'hgvsg',
            'hgvsc',
            'hgvsp',
            'hgnc_id',
            'amino_acids',
            'impact',
            'consequence_terms',
            'codons',
            'biotype',
            'polyphen_score',
            'strand',
                'sift_score']:
            output_item(transcript, id, indent, highlight)
        output_kv('cds_position', format_position(
            transcript, 'cds'), indent, highlight)
        output_kv('cdna_position', format_position(
            transcript, 'cdna'), indent, highlight)
        output_kv('protein_position', format_protein_position(
            transcript), indent, highlight)
        output_kv('transcript_name', transcript_name(
            transcript), indent, highlight)

    def format_allele_string(allele_string):
        return "%s > %s" % tuple(allele_string.split('/'))

    def format_protein_position(transcript):
        if 'protein_start' in transcript and 'protein_end' in transcript:
            return '%s-%s' % (transcript['protein_start'], transcript['protein_end']) if transcript['protein_start'] != transcript['protein_end'] else transcript['protein_start']
        else:
            return 'N/A'

    def format_position(transcript, prefix):
        if prefix+'_start' in transcript and prefix+'_end' in transcript:
            if transcript[prefix+'_start'] == transcript[prefix+'_end']:
                return "%s" % (transcript[prefix+'_start'])
            else:
                return "%s-%s" % (transcript[prefix+'_start'], transcript[prefix+'_end'])
        else:
            return 'N/A'

    def format_change_position(variant):
        if 'start' in variant and 'end' in variant:
            return '%s-%s' % (variant['start'], variant['end']) if variant['start'] != variant['end'] else variant['start']
        else:
            return 'N/A'

    def transcript_name(transcript):
        if 'amino_acids' in transcript:
            amino_from, amino_to = transcript['amino_acids'].split('/')
            position = format_protein_position(transcript)
            return '%s%s%s' % (amino_from, position, amino_to)
        else:
            return 'N/A'

    def enrich_transcript(transcript, species=species, GRCh37=GRCh37):
        if 'hgvsc' in transcript:
            from bioinfo_toolset.modules.variant_recoder import recode
            records = recode(transcript['hgvsc'], species, GRCh37)
            transcript['hgvsg'] = []
            for rec in records:
                for allele in rec.values():
                    for hgvsg in allele['hgvsg']:
                        transcript['hgvsg'].append(hgvsg)
        return transcript

    try:
        if input_type == 'vcf':
            offline_vep = OfflineVep()
            offline_vep.populate_cache()
            results = offline_vep.evaluate(
                [input],
                GRCh37=GRCh37
            )
        elif input_type == 'vcf_file':
            offline_vep = OfflineVep()
            offline_vep.populate_cache()
            lines = []
            with open(input, 'r') as inf:
                for line in inf.readlines():
                    if not line.startswith('#'):
                        lines.append(line.strip())
            results = offline_vep.evaluate(
                lines,
                GRCh37=GRCh37
            )
        else:
            results = vep(input, species=species,
                          input_type=input_type, GRCh37=GRCh37, refseq=refsec_mode)

        for result in results:
            print(
                colored(f"Result for {result['input']}", 'blue', attrs={'bold'}))
            print(
                colored(f"{result['assembly_name']}: {result['id']} -> {result['seq_region_name']}:{format_change_position(result)} ({format_allele_string(result['allele_string'])})", 'green', attrs=['bold']))
            if _liftover:
                lo_postion = {}
                t_chr, lo_postion['start'] = liftover(
                    'hg19' if GRCh37 else 'hg38', 'hg38' if GRCh37 else 'hg19', result['seq_region_name'], result['start'])
                t_chr, lo_postion['end'] = liftover(
                    'hg19' if GRCh37 else 'hg38', 'hg38' if GRCh37 else 'hg19', result['seq_region_name'], result['end'])
                print(
                    colored(f"{'GRCh38' if GRCh37 else 'GRCh37'}: {result['id']} -> {t_chr}:{format_change_position(lo_postion)} ({format_allele_string(result['allele_string'])})", 'green', attrs=['bold']))
            output_variant(result)

            if 'colocated_variants' in result:
                print(
                    colored(f"Colocated variants ({len(result['colocated_variants'])}):", 'cyan'))
                for colocated_variant in result['colocated_variants']:
                    output_variant(colocated_variant, indent=1)
                    output_splitter(indent=1)

            if 'transcript_consequences' in result:
                print(
                    colored(f"Transcript consequences ({len(result['transcript_consequences'])}{' / only canonical shown' if not all_transcripts else ''}):", 'cyan'))
                for transcript_consequence in result['transcript_consequences']:
                    if all_transcripts or (not all_transcripts and 'canonical' in transcript_consequence and transcript_consequence['canonical'] == 1):
                        if enrich_transcripts:
                            transcript_consequence = enrich_transcript(
                                transcript_consequence)
                        output_transcript(transcript_consequence, indent=1)
                        output_splitter(indent=1)

    except Exception as ex:
        log.error(ex)
        traceback.print_exc(file=stdout)


cli.add_command(vep)

if __name__ == '__main__':
    cli()
