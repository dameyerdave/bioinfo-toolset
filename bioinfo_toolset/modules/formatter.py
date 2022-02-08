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
    if 'amino_acids' in transcript:
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
    else:
        return suggestion
