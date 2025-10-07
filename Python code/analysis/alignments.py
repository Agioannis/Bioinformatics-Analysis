from Bio.Align import PairwiseAligner

def sequence_alignment(seq1, seq2, method="nw"):
    seq1 = str(seq1)
    seq2 = str(seq2)

    aligner = PairwiseAligner()

    # Set scoring scheme
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    if method == "nw":  # Needleman-Wunsch (global)
        aligner.mode = "global"
    elif method == "sw":  # Smith-Waterman (local)
        aligner.mode = "local"
    else:
        raise ValueError("Invalid method. Use 'nw' or 'sw'.")

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return None

    aln = alignments[0]  # best alignment
    aligned_seq1 = aln.aligned[0]
    aligned_seq2 = aln.aligned[1]

    # Convert alignment into readable strings
    alignment_str = str(aln)

    # Identity calculation (rough: compare aligned sequences)
    aligned_str1, aligned_str2 = alignment_str.split("\n")[:2]
    matches = sum(a == b for a, b in zip(aligned_str1, aligned_str2) if a != "-" and b != "-")
    identity = (matches / min(len(seq1), len(seq2))) * 100

    # Count gaps
    gaps = aligned_str1.count("-") + aligned_str2.count("-")

    return {
        "score": aln.score,
        "identity": identity,
        "gaps": gaps,
        "alignment": alignment_str
    }
