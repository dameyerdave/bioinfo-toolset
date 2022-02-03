def format_allele_string(allele_string):
    if '/' in allele_string:
        return "%s > %s" % tuple(allele_string.split('/'))
    else:
        return allele_string


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
