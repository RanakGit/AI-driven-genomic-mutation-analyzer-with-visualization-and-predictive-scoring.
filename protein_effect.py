from Bio.Seq import Seq
from Bio.Data import CodonTable


def translate_sequence(dna_seq, frame=0):
    seq = Seq(dna_seq[frame:])
    return str(seq.translate(to_stop=True))


def variant_protein_impact(ref_seq, variant):
    """Assess variant impact based on codon/protein change heuristics."""
    # variant has keys: type, position_ref, ref, alt
    ptype = variant['type']
    pos = variant['position_ref']

    if ptype == 'SNP':
        codon_idx = (pos - 1) // 3
        codon_start = codon_idx * 3
        if codon_start + 3 > len(ref_seq):
            return 'Unknown'

        ref_codon = ref_seq[codon_start:codon_start + 3]
        alt_seq = list(ref_seq)
        alt_seq[pos - 1] = variant['alt']
        alt_seq = ''.join(alt_seq)
        alt_codon = alt_seq[codon_start:codon_start + 3]

        if len(ref_codon) != 3 or len(alt_codon) != 3:
            return 'Unknown'

        try:
            ref_aa = Seq(ref_codon).translate()
            alt_aa = Seq(alt_codon).translate()
        except Exception:
            return 'Unknown'

        if ref_aa == alt_aa:
            return 'Synonymous'
        if ref_aa == '*' or alt_aa == '*':
            return 'Nonsense'
        return 'Missense'

    if ptype in ('INS', 'DEL'):
        return 'Frameshift' if (len(variant['alt']) - len(variant['ref'])) % 3 != 0 else 'In-frame'

    return 'Unknown'
