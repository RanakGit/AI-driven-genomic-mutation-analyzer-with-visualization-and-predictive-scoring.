try:
    from Bio.Align import PairwiseAligner
    ALIGNER_AVAILABLE = True
except ImportError:
    from Bio import pairwise2
    ALIGNER_AVAILABLE = False


def align_sequences(seq1: str, seq2: str):
    # Fast direct diff if sequences are the same length with no indels.
    if len(seq1) == len(seq2):
        score = sum(2 if a == b else -1 for a, b in zip(seq1, seq2))
        return seq1, seq2, score, 0, len(seq1)

    if ALIGNER_AVAILABLE:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        alignment = aligner.align(seq1, seq2)[0]
        seqA, seqB = alignment.sequences
        return str(seqA), str(seqB), alignment.score, 0, len(seqA)

    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
    if len(alignments) == 0:
        raise ValueError("No alignment found")
    return alignments[0]


def detect_variants(seq1: str, seq2: str):
    # seq1: normal/reference, seq2: query/mutant
    aln1, aln2, score, start, end = align_sequences(seq1, seq2)
    variants = []
    pos_ref = 0
    pos_query = 0

    for a, b in zip(aln1, aln2):
        if a != '-' and b != '-':
            pos_ref += 1
            pos_query += 1
            if a != b:
                variants.append({
                    'type': 'SNP',
                    'position_ref': pos_ref,
                    'position_query': pos_query,
                    'ref': a,
                    'alt': b,
                })
        elif a == '-' and b != '-':
            pos_query += 1
            variants.append({
                'type': 'INS',
                'position_ref': pos_ref,
                'position_query': pos_query,
                'ref': '-',
                'alt': b,
            })
        elif a != '-' and b == '-':
            pos_ref += 1
            variants.append({
                'type': 'DEL',
                'position_ref': pos_ref,
                'position_query': pos_query,
                'ref': a,
                'alt': '-',
            })

    return {
        'alignment_ref': aln1,
        'alignment_query': aln2,
        'score': score,
        'variants': variants
    }


def summarize_variants(variants):
    counts = {'SNP': 0, 'INS': 0, 'DEL': 0}
    for v in variants:
        counts[v['type']] = counts.get(v['type'], 0) + 1
    return counts
