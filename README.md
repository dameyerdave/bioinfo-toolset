# BioInformatics Tools

## Installation

```bash
python3 -m pip install bioinfo-toolset
```

## Usage

```bash
bit --help
```

## Example

```bash
bit vep -l ENST00000257430.9:c.3927_3931del
```

```bash
Result for ENST00000257430.9:c.3927_3931del
GRCh38: ENST00000257430.9:c.3927_3931del -> 5:112839521-112839525 (AAAGA > -)
GRCh37: ENST00000257430.9:c.3927_3931del -> 5:112175218-112175222 (AAAGA > -)
 assembly_name         : GRCh38
 id                    : ENST00000257430.9:c.3927_3931del
 seq_region_name       : 5
 start                 : 112839521
 end                   : 112839525
 allele_string         : AAAGA/-
 strand                : 1
 most_severe_consequence : frameshift_variant
Colocated variants (1):
    id                    : COSV57321812
    seq_region_name       : 5
    start                 : 112839521
    end                   : 112839525
    allele_string         : COSMIC_MUTATION
    strand                : 1
    somatic               : 1
    ---
Transcript consequences (8):
    transcript_id         : ENST00000257430
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    amino_acids           : EKI/DX
    consequence_terms     : ['frameshift_variant']
    protein_position      : 1309-1311
    transcript_name       : EKI1309-1311DX
    ---
    transcript_id         : ENST00000502371
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    consequence_terms     : ['3_prime_UTR_variant', 'NMD_transcript_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
    transcript_id         : ENST00000504915
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    consequence_terms     : ['downstream_gene_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
    transcript_id         : ENST00000507379
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    consequence_terms     : ['downstream_gene_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
    transcript_id         : ENST00000508376
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    amino_acids           : EKI/DX
    consequence_terms     : ['frameshift_variant']
    protein_position      : 1309-1311
    transcript_name       : EKI1309-1311DX
    ---
    transcript_id         : ENST00000508624
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    consequence_terms     : ['3_prime_UTR_variant', 'NMD_transcript_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
    transcript_id         : ENST00000512211
    gene_symbol           : APC
    gene_id               : ENSG00000134982
    hgnc_id               : HGNC:583
    consequence_terms     : ['downstream_gene_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
    transcript_id         : ENST00000520401
    gene_id               : ENSG00000258864
    consequence_terms     : ['intron_variant', 'NMD_transcript_variant']
    protein_position      : N/A
    transcript_name       : N/A
    ---
```