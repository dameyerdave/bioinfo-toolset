#!/usr/bin/env python

import click
from click import Choice
from friendlylog import colored_logger as log
from termcolor import colored
import traceback
from sys import stdout


@click.group()
def cli():
    pass


@click.command()
@click.option('--from', '-f', '_from', default='hg19', type=str, help='Convert from this genome assembly')
@click.option('--to', '-t', default='hg38', type=str, help='Convert to this genome assembly')
@click.argument('chromosome', type=int)
@click.argument('position', type=int)
def liftover(_from, to, chromosome, position):
    from bioinfo_toolset.liftover import liftover
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


@click.command(context_settings={"ignore_unknown_options": True}, help="The query string can be a hgvs, id or region (e.g. hgvs: 9:g.22125504G>C, region: 3:178928079/C)")
@click.option('--GRCh37', '--grch37', '--old', '-o', 'GRCh37', is_flag=True, help="Use the GRCh rest api for the query")
@click.option('--liftover', '-lift', '-l', 'liftover', is_flag=True, help="Liftover the variant GRCh38/GRCh37, depending on the used assembly")
@click.option('--species', '-s', default='human', type=str, help='The species to query VEP')
@click.option('--input-type', '-t', 'input_type', default='hgvs', type=Choice(['hgvs', 'id', 'region']), help='The input type you want to query')
@click.argument('input', type=str)
def vep(species, input_type, input, GRCh37, liftover):
    from bioinfo_toolset.vep import vep
    from bioinfo_toolset.liftover import liftover

    def output(text, indent=0):
        print(f'%-{indent * 4}s{text}' % (' '))

    def output_kv(id, value, indent=0, highlight=False):
        attrs = []
        if highlight:
            attrs = ['bold']
        print(f'%-{indent * 4}s%-{34 if highlight else 30}s : %s' %
              (' ', colored(id, 'cyan', attrs=attrs), colored(value, 'white', attrs=attrs)))

    def output_item(_dict, id, indent=0, highlight=False):
        if id in _dict:
            output_kv(id, _dict[id], indent, highlight)

    def output_splitter(indent=0):
        output('---', indent)

    def output_variant(variant, indent=0):
        for id in ['assembly_name', 'id', 'seq_region_name', 'start', 'end', 'allele_string', 'strand', 'somatic', 'most_severe_consequence']:
            output_item(variant, id, indent)

    def output_transcript(transcript, indent=0):
        highlight = False
        if 'canonical' in transcript and transcript['canonical'] == 1:
            highlight = True

        for id in ['transcript_id', 'gene_symbol', 'gene_id', 'hgnc_id', 'amino_acids', 'protein_start' 'impact', 'consequence_terms']:
            output_item(transcript, id, indent, highlight)
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

    try:
        results = vep(input, species=species,
                      input_type=input_type, GRCh37=GRCh37)
        for result in results:
            print(
                colored(f"Result for {result['input']}", 'blue', attrs={'bold'}))
            print(
                colored(f"{result['assembly_name']}: {result['id']} -> {result['seq_region_name']}:{format_change_position(result)} ({format_allele_string(result['allele_string'])})", 'green', attrs=['bold']))
            if liftover:
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
                    colored(f"Transcript consequences ({len(result['transcript_consequences'])}):", 'cyan'))
                for transcript_consequence in result['transcript_consequences']:
                    output_transcript(transcript_consequence, indent=1)
                    output_splitter(indent=1)

    except Exception as ex:
        log.error(ex)
        traceback.print_exc(file=stdout)


cli.add_command(vep)

if __name__ == '__main__':
    cli()
