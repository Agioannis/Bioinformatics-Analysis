import logging
import re
from typing import List, Dict
from Bio.Seq import Seq

def _bases_match(base1: str, base2: str) -> bool:
    # Ελέγχει αν δύο βάσεις ταιριάζουν (N ταιριάζει με οποιαδήποτε βάση)
    if base1 == 'N' or base2 == 'N':
        return True  # Το N ταιριάζει με οποιαδήποτε βάση
    return base1 == base2

def _patterns_match(pattern1: str, pattern2: str) -> bool:
    # Ελέγχει αν δύο μοτίβα νουκλεοτιδίων ταιριάζουν θέση προς θέση
    if len(pattern1) != len(pattern2):
        return False
    return all(_bases_match(b1, b2) for b1, b2 in zip(pattern1, pattern2))

def find_tandem_repeats(seq: Seq, min_period=2, max_period=10, min_repeats=3, max_results=50) -> List[Dict]:
    # Βρίσκει τα ενωμένα επαναλαμβανόμενα μοτίβα (tandem repeats) σε μια ακολουθία DNA
    try:
        sequence = str(seq).upper()
        repeats = []

        # Αγνοεί πολύ μικρές ακολουθίες
        if len(sequence) < min_period * min_repeats:
            return []

        # Αναζητά επαναλήψεις για κάθε φάσμα περιόδου
        for period in range(min_period, min(max_period + 1, len(sequence) // min_repeats + 1)):
            # μοτίβο με παράθυρο sliding
            for start in range(len(sequence) - period * min_repeats + 1):
                pattern = sequence[start:start + period]

                # Μετρά διαδοχικές επαναλήψεις
                repeat_count = 1
                pos = start + period

                while pos + period <= len(sequence) and _patterns_match(sequence[pos:pos + period], pattern):
                    repeat_count += 1
                    pos += period

                # Ελέγχει αν πληροί τα κριτήρια ελάχιστων επαναλήψεων
                if repeat_count >= min_repeats:
                    end_pos = start + (repeat_count * period)

                    # Ελέγχει για επικαλύψεις με υπάρχουσες επαναλήψεις
                    overlap = False
                    for existing in repeats:
                        if not (end_pos <= existing['start'] or start >= existing['end']):
                            # Αν η νέα επανάληψη είναι μεγαλύτερη, αφαιρεί την παλιά
                            if (end_pos - start) > (existing['end'] - existing['start']):
                                repeats.remove(existing)
                                break
                            else:
                                overlap = True
                                break

                    if not overlap:
                        repeat_info = {
                            'pattern': pattern,
                            'start': start + 1,
                            'end': end_pos,
                            'repeat_count': repeat_count,
                            'total_length': end_pos - start,
                            'period': period,
                            'sequence': sequence[start:end_pos],
                            'purity': 1.0,  # Τέλειες επαναλήψεις
                            'type': _classify_repeat_type(pattern, period)
                        }
                        repeats.append(repeat_info)

        # Ταξινόμηση κατά συνολικό μήκος (φθίνουσα) και μετά κατά θέση έναρξης
        repeats.sort(key=lambda x: (-x['total_length'], x['start']))

        return repeats[:max_results]

    except Exception as e:
        logging.warning(f"Ανάλυση επαναλήψεων απέτυχε: {e}")
        return []

def find_imperfect_repeats(seq: Seq, min_period=2, max_period=6, min_repeats=3, max_mismatches=1) -> List[Dict]:
    # Βρίσκει μη τέλειες επαναλήψεις με επιτρεπτές αποκλίσεις (mismatches)
    try:
        sequence = str(seq).upper()
        imperfect_repeats = []

        if len(sequence) < min_period * min_repeats:
            return []

        for period in range(min_period, min(max_period + 1, len(sequence) // min_repeats + 1)):
            for start in range(len(sequence) - period * min_repeats + 1):
                pattern = sequence[start:start + period]

                repeat_count = 1
                total_mismatches = 0
                pos = start + period

                while pos + period <= len(sequence):
                    current_unit = sequence[pos:pos + period]
                    mismatches = sum(1 for i in range(period) if not _bases_match(pattern[i], current_unit[i]))

                    if mismatches <= max_mismatches:
                        repeat_count += 1
                        total_mismatches += mismatches
                        pos += period
                    else:
                        break

                if repeat_count >= min_repeats:
                    end_pos = start + (repeat_count * period)
                    purity = 1.0 - (total_mismatches / (repeat_count * period))

                    repeat_info = {
                        'pattern': pattern,
                        'start': start + 1,
                        'end': end_pos,
                        'repeat_count': repeat_count,
                        'total_length': end_pos - start,
                        'period': period,
                        'sequence': sequence[start:end_pos],
                        'purity': purity,
                        'mismatches': total_mismatches,
                        'type': _classify_repeat_type(pattern, period)
                    }
                    imperfect_repeats.append(repeat_info)

        # Αφαιρεί επικαλύψεις και ταξινομεί
        imperfect_repeats = _remove_overlapping_repeats(imperfect_repeats)
        imperfect_repeats.sort(key=lambda x: (-x['total_length'], x['start']))

        return imperfect_repeats[:20]

    except Exception as e:
        logging.warning(f"Ανάλυση μη τέλειων επαναλήψεων απέτυχε: {e}")
        return []

def find_microsatellites(seq: Seq, min_repeats=4) -> List[Dict]:
    try:
        sequence = str(seq).upper()
        microsatellites = []

        # Κοινά μοτίβα μικροδορυφόρων ανά μήκος μοτίβου
        patterns = {
            1: ['A', 'T', 'G', 'C'],
            2: ['AT', 'TA', 'GC', 'CG', 'AG', 'GA', 'CT', 'TC', 'AC', 'CA', 'TG', 'GT'],
            3: ['ATG', 'TGA', 'GAT', 'CAT', 'ATC', 'TCA', 'GCA', 'CAG', 'AGC', 'GCT', 'CTG', 'TGC'],
            4: ['ATGC', 'TGCA', 'GCAT', 'CATG', 'AAAG', 'AAGA', 'AGAA', 'GAAA', 'TTTC', 'TTCT', 'TCTT', 'CTTT'],
            5: ['ATGCA', 'TGCAT', 'GCATG', 'CATGA', 'ATGAC'],
            6: ['ATGCAT', 'TGCATG', 'GCATGA', 'CATGAT']
        }

        for period in range(1, 7):
            for pattern in patterns[period]:
                regex_pattern = f'({re.escape(pattern)}){{{min_repeats},}}'
                matches = re.finditer(regex_pattern, sequence, re.IGNORECASE)

                for match in matches:
                    start_pos = match.start()
                    end_pos = match.end()
                    matched_seq = match.group()
                    repeat_count = len(matched_seq) // period

                    microsatellite_info = {
                        'pattern': pattern,
                        'start': start_pos + 1,
                        'end': end_pos,
                        'repeat_count': repeat_count,
                        'total_length': end_pos - start_pos,
                        'period': period,
                        'sequence': matched_seq,
                        'purity': 1.0,
                        'type': f"Microsatellite ({period}-mer)"
                    }
                    microsatellites.append(microsatellite_info)

        # Αφαιρεί επικαλύψεις και ταξινομεί
        microsatellites = _remove_overlapping_repeats(microsatellites)
        microsatellites.sort(key=lambda x: (-x['total_length'], x['start']))

        return microsatellites[:30]

    except Exception as e:
        logging.warning(f"Ανάλυση μικροδορυφόρων απέτυχε: {e}")
        return []

def _classify_repeat_type(pattern: str, period: int) -> str:
    # Κατηγοριοποιεί τον τύπο της επανάληψης βάσει του μοτίβου και της περιόδου
    if period == 1:
        return f"Ομοπολυμερές ({pattern})"
    elif period == 2:
        return f"Δι-νουκλεοτιδικό ({pattern})"
    elif period == 3:
        return f"Τρι-νουκλεοτιδικό ({pattern})"
    elif period == 4:
        return f"Τετρα-νουκλεοτιδικό ({pattern})"
    elif period == 5:
        return f"Πεντα-νουκλεοτιδικό ({pattern})"
    elif period == 6:
        return f"Εξα-νουκλεοτιδικό ({pattern})"
    else:
        return f"Σύνθετο ({period}-mer)"

def _remove_overlapping_repeats(repeats: List[Dict]) -> List[Dict]:
    # Αφαιρεί επικαλυπτόμενες επαναλήψεις, διατηρώντας τις μεγαλύτερες
    if not repeats:
        return []

    sorted_repeats = sorted(repeats, key=lambda x: -x['total_length'])
    non_overlapping = []

    for repeat in sorted_repeats:
        overlap = False
        for existing in non_overlapping:
            if not (repeat['end'] <= existing['start'] or repeat['start'] >= existing['end']):
                overlap = True
                break

        if not overlap:
            non_overlapping.append(repeat)

    return non_overlapping

def analyze_repeat_distribution(seq: Seq) -> Dict:
    # Αναλύει τη συνολική κατανομή των επαναλήψεων σε μια ακολουθία
    try:
        all_repeats = find_tandem_repeats(seq, max_results=100)
        microsatellites = find_microsatellites(seq)

        if not all_repeats:
            return {
                'total_repeats': 0,
                'repeat_coverage': 0.0,
                'average_repeat_length': 0.0,
                'repeat_density': 0.0,
                'most_common_period': 0,
                'repeat_types': {}
            }

        seq_length = len(str(seq))
        total_repeat_length = sum(repeat['total_length'] for repeat in all_repeats)

        repeat_types = {}
        period_counts = {}

        for repeat in all_repeats:
            repeat_type = repeat['type']
            period = repeat['period']

            if repeat_type not in repeat_types:
                repeat_types[repeat_type] = 0
            repeat_types[repeat_type] += 1

            if period not in period_counts:
                period_counts[period] = 0
            period_counts[period] += 1

        most_common_period = max(period_counts.keys(), key=lambda k: period_counts[k]) if period_counts else 0

        return {
            'total_repeats': len(all_repeats),
            'repeat_coverage': (total_repeat_length / seq_length) * 100,
            'average_repeat_length': total_repeat_length / len(all_repeats),
            'repeat_density': len(all_repeats) / (seq_length / 1000),  # ανά kb
            'most_common_period': most_common_period,
            'repeat_types': repeat_types,
            'microsatellites_found': len(microsatellites)
        }

    except Exception as e:
        logging.warning(f"Ανάλυση κατανομής επαναλήψεων απέτυχε: {e}")
        return {
            'total_repeats': 0,
            'repeat_coverage': 0.0,
            'average_repeat_length': 0.0,
            'repeat_density': 0.0,
            'most_common_period': 0,
            'repeat_types': {},
            'microsatellites_found': 0
        }
