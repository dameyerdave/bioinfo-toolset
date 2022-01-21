_complement_allele_lookup = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}


def complement_allele_lookup(reference_allele):
    complement_allele = ''
    for a in reference_allele:
        complement_allele += _complement_allele_lookup[a]
    return complement_allele
