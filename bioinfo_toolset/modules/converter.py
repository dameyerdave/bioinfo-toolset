from Bio.PDB.Polypeptide import three_to_one as biopython_three_to_one, standard_aa_names
import re

def three_to_one(three_letter_string: str):
    converted = three_letter_string
    for aa_name in standard_aa_names:
        converted = re.sub(aa_name, biopython_three_to_one(aa_name), converted, flags=re.I)
    return converted

