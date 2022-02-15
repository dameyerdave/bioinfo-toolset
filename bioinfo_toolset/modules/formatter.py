from bioinfo_toolset.modules.converter import three_to_one
import re
from friendlylog import colored_logger as log


def format_allele_string(allele_string):
    if '/' in allele_string:
        return "%s > %s" % tuple(allele_string.split('/'))
    else:
        return allele_string


def format_protein_position(transcript, delim='-'):
    if 'protein_start' in transcript and 'protein_end' in transcript:
        return '%s%s%s' % (transcript['protein_start'], delim, transcript['protein_end']) if transcript['protein_start'] != transcript['protein_end'] else transcript['protein_start']
    else:
        return None


def format_position(transcript, prefix):
    if prefix+'_start' in transcript and prefix+'_end' in transcript:
        if transcript[prefix+'_start'] == transcript[prefix+'_end']:
            return "%s" % (transcript[prefix+'_start'])
        else:
            return "%s-%s" % (transcript[prefix+'_start'], transcript[prefix+'_end'])
    else:
        return None


def format_change_position(variant):
    if 'start' in variant and 'end' in variant:
        return '%s-%s' % (variant['start'], variant['end']) if variant['start'] != variant['end'] else variant['start']
    else:
        return None


variant_allele_lookup = {
    'deletion': 'del',
    'duplication': 'dup'
}


def transcript_name(transcript, suggestion=None):
    """Returns the transcript name of this transcript in one char nomenclature"""
    try:
        # First we try to three_to_one the hgvsp
        if 'hgvsp' in transcript:
            return re.sub(r'p\.', '', three_to_one(transcript['hgvsp'].split(':')[1]))
            # If the hgvsp is not given we try to combine amino accids to
            # create a name

        elif 'amino_acids' in transcript:
            if '/' in transcript['amino_acids']:
                amino_from, amino_to = transcript['amino_acids'].split('/')
            else:
                amino_from = ''
                if transcript['variant_allele'] in variant_allele_lookup:
                    amino_to = variant_allele_lookup[transcript['variant_allele']]
                else:
                    # this is only to debug purposes to fill the variant_allele_lookup appropriately
                    amino_to = transcript['variant_allele']
            position = format_protein_position(transcript, delim='_')
            return '%s%s%s' % (amino_from, position, amino_to)
        # Splice site variants
        elif 'consequence_terms' in transcript and 'hgvsc' in transcript:
            if 'splice_acceptor_variant' in transcript['consequence_terms']:
                return "splice site {}".format(re.sub(r'[^:]+:[a-z]\.(.+)$', r'\1', transcript['hgvsc']))
    except:
        log.warning(
            f"Unable to build transcript name for {transcript['hgvsg'] if 'hgvsg' in transcript else transcript['transcript_id']}.")

    # if no name can be found we return the suggestion or None if no suggestion is given
    return suggestion
