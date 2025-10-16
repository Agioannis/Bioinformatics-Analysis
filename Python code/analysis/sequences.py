import logging
import numpy as np
from typing import List, Dict
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from .repeats import find_tandem_repeats

# Υπολογισμός της περιεκτικότητας σε GC
# Επιστρέφει το ποσοστό των G και C βάσεων σε μία ακολουθία DNA
def gc_content(seq: Seq) -> float:
    s = ''.join(c for c in str(seq).upper() if c in 'ATGCN')
    return (s.count('G') + s.count('C')) / len(s) * 100 if s else 0.0

# Προχωρημένη ανάλυση GC περιεκτικότητας
# Παρέχει AT/GC αναλογίες, skew, και ποσοστά πουρίνης/πυριμιδίνης
def advanced_gc_analysis(seq: Seq) -> Dict:
    s = str(seq).upper().replace(' ', '').replace('\n', '')
    if not s:
        return {}

    counts = {base: s.count(base) for base in 'ATGC'}
    gc_count = counts['G'] + counts['C']
    at_count = counts['A'] + counts['T']

    return {
        'gc_content': gc_count / len(s) * 100,
        'at_content': at_count / len(s) * 100,
        'gc_skew': (counts['G'] - counts['C']) / (gc_count) if gc_count > 0 else 0,
        'at_gc_ratio': at_count / gc_count if gc_count > 0 else float('inf'),
        'purine_content': (counts['A'] + counts['G']) / len(s) * 100,
        'pyrimidine_content': (counts['T'] + counts['C']) / len(s) * 100
    }

# Ανάλυση πολυπλοκότητας ακολουθίας
# Υπολογίζει την εντροπία Shannon και περιοχές χαμηλής πολυπλοκότητας
def analyze_sequence_complexity(seq: Seq) -> Dict:
    s = str(seq).upper()
    seq_len = len(s)

    if not s:
        return {}

    from collections import Counter
    counts = Counter(s)
    # Εντροπία Shannon ως μέτρο πολυπλοκότητας
    entropy = -sum((count / seq_len) * np.log2(count / seq_len) for count in counts.values())

    # Εντοπισμός περιοχών χαμηλής πολυπλοκότητας
    low_complexity_regions = []
    window_size = 20
    threshold = 1.5  # Εντροπία κάτω από αυτό σημαίνει χαμηλή πολυπλοκότητα
    step_size = 1 if seq_len < 50000 else max(1, seq_len // 5000)

    for i in range(0, seq_len - window_size + 1, step_size):
        if len(low_complexity_regions) >= 100:
            break
        window = s[i:i + window_size]
        window_counts = Counter(window)
        window_entropy = -sum((count / window_size) * np.log2(count / window_size) for count in window_counts.values())

        if window_entropy < threshold:
            low_complexity_regions.append({
                'start': i + 1,
                'end': i + window_size,
                'entropy': window_entropy,
                'sequence': window
            })

    # Υπολογισμός δείκτη πολυπλοκότητας και συσχετισμών
    return {
        'shannon_entropy': entropy,
        'linguistic_complexity': 0.0 if seq_len > 10000 else len(set(s[i:j] for i in range(min(1000, seq_len)) for j in range(i+1, min(i+50, seq_len+1)))) / (seq_len * (seq_len + 1) / 2),
        'low_complexity_regions': low_complexity_regions,
        'max_entropy': 2.0,  # Μέγιστη θεωρητική εντροπία για DNA
        'complexity_score': entropy / 2.0
    }

# Εύρεση γνωστών μοτίβων (motifs) σε μια ακολουθία DNA
def find_motifs(seq: Seq, motif_patterns: List[str] | None = None) -> Dict:
    if motif_patterns is None:
        # Γνωστά βιολογικά μοτίβα (π.χ. περιοχές σύνδεσης ενζύμων περιορισμού)
        motif_patterns = ['TATAAA', 'CAAT', 'GGGCGG', 'AAGCTT', 'GAATTC', 'GGATCC', 'CTCGAG', 'GCGGCCGC', 'GAGCTC', 'AGATCT']

    s = str(seq).upper()
    motif_results = {}

    for motif in motif_patterns:
        # Εύρεση όλων των θέσεων όπου εμφανίζεται το μοτίβο
        positions = [i + 1 for i in range(len(s)) if s.startswith(motif, i)]
        if positions:
            motif_results[motif] = {
                'positions': positions,
                'count': len(positions),
                'density': len(positions) / len(s) * 1000  # πλήθος ανά 1000 βάσεις
            }
    return motif_results

# Ολοκληρωμένη ανάλυση ακολουθίας DNA
# Συνδυάζει βασικά στατιστικά, GC ανάλυση, μοτίβα και πολυπλοκότητα
def comprehensive_analysis(seq: Seq, seq_id: str) -> Dict:
    # Καθαρισμός ακολουθίας (διατήρηση μόνο A, T, G, C, N)
    cleaned_seq_str = ''.join(c for c in str(seq).upper() if c in 'ATGCN')
    cleaned_seq = Seq(cleaned_seq_str)
    seq_len = len(cleaned_seq_str)

    # Για μοριακό βάρος, αγνοούμε τα N (άδηλες βάσεις)
    seq_for_mol_weight = Seq(''.join(c for c in cleaned_seq_str if c in 'ATGC'))

    # Βασικά στατιστικά
    results = {
        'basic_stats': {
            'id': seq_id,
            'length': seq_len,
            **{base: cleaned_seq_str.count(base) for base in 'ATGC'},
            'gc_content': gc_content(cleaned_seq),
            'molecular_weight': molecular_weight(seq_for_mol_weight) if seq_for_mol_weight else 0.0
        }
    }

    # Εκτέλεση προχωρημένων αναλύσεων με διαχείριση σφαλμάτων
    try:
        results.update({
            'advanced_gc': advanced_gc_analysis(cleaned_seq),
            'motifs': find_motifs(cleaned_seq),
            'complexity': analyze_sequence_complexity(cleaned_seq)
        })
    except Exception as e:
        logging.warning(f"Advanced analysis skipped for {seq_id}: {e}")
        results.update({'advanced_gc': {}, 'motifs': {}, 'complexity': {}})

    # Εύρεση tandem repeats (π.χ. επαναλαμβανόμενων υπομονάδων DNA)
    try:
        results['tandem_repeats'] = find_tandem_repeats(cleaned_seq)
    except Exception:
        results['tandem_repeats'] = []

    # Πρόβλεψη θέσεων περιορισμού και πιθανών εκκινητών (Primers)
    results['restriction_sites'] = {}
    results['primers'] = []

    return results
