
import logging
import re
from typing import List, Dict
from Bio.Seq import Seq

def _bases_match(base1: str, base2: str) -> bool:
    """
    Check if two DNA bases match, treating N as a wildcard

    Parameters:
    - base1, base2: Single nucleotide characters

    Returns:
    True if bases match (including when either is N), False otherwise
    """
    if base1 == 'N' or base2 == 'N':
        return True  # N matches any base
    return base1 == base2

def _patterns_match(pattern1: str, pattern2: str) -> bool:
    """
    Check if two DNA patterns match, treating N as wildcards

    Parameters:
    - pattern1, pattern2: DNA sequence strings of same length

    Returns:
    True if patterns match (with N as wildcards), False otherwise
    """
    if len(pattern1) != len(pattern2):
        return False
    return all(_bases_match(b1, b2) for b1, b2 in zip(pattern1, pattern2))

def find_tandem_repeats(seq: Seq, min_period=2, max_period=10, min_repeats=3, max_results=50) -> List[Dict]:
    """
    Find tandem repeats in a DNA sequence

    Parameters:
    - seq: BioPython Seq object
    - min_period: minimum repeat unit length (default 2)
    - max_period: maximum repeat unit length (default 10)
    - min_repeats: minimum number of consecutive repeats (default 3)
    - max_results: maximum number of results to return (default 50)

    Returns:
    List of dictionaries containing repeat information
    """
    try:
        sequence = str(seq).upper()
        repeats = []

        # Skip very short sequences
        if len(sequence) < min_period * min_repeats:
            return []

        # Find repeats for each period length
        for period in range(min_period, min(max_period + 1, len(sequence) // min_repeats + 1)):
            # Sliding window approach
            for start in range(len(sequence) - period * min_repeats + 1):
                pattern = sequence[start:start + period]

                # Handle patterns with ambiguous bases (N)
                # Instead of skipping, we'll treat N as a wildcard that can match any base

                # Count consecutive repeats
                repeat_count = 1
                pos = start + period

                while pos + period <= len(sequence) and _patterns_match(sequence[pos:pos + period], pattern):
                    repeat_count += 1
                    pos += period

                # Check if meets minimum repeat criteria
                if repeat_count >= min_repeats:
                    end_pos = start + (repeat_count * period)

                    # Check for overlaps with existing repeats
                    overlap = False
                    for existing in repeats:
                        if not (end_pos <= existing['start'] or start >= existing['end']):
                            # If new repeat is longer, remove the existing one
                            if (end_pos - start) > (existing['end'] - existing['start']):
                                repeats.remove(existing)
                                break
                            else:
                                overlap = True
                                break

                    if not overlap:
                        repeat_info = {
                            'pattern': pattern,
                            'start': start + 1,  # 1-based indexing
                            'end': end_pos,
                            'repeat_count': repeat_count,
                            'total_length': end_pos - start,
                            'period': period,
                            'sequence': sequence[start:end_pos],
                            'purity': 1.0,  # Perfect repeats have 100% purity
                            'type': _classify_repeat_type(pattern, period)
                        }
                        repeats.append(repeat_info)

        # Sort by total length (descending) and then by start position
        repeats.sort(key=lambda x: (-x['total_length'], x['start']))

        # Return top results
        return repeats[:max_results]

    except Exception as e:
        logging.warning(f"Tandem repeat analysis failed: {e}")
        return []

def find_imperfect_repeats(seq: Seq, min_period=2, max_period=6, min_repeats=3, max_mismatches=1) -> List[Dict]:
    """
    Find imperfect tandem repeats allowing for mismatches

    Parameters:
    - seq: BioPython Seq object
    - min_period: minimum repeat unit length
    - max_period: maximum repeat unit length
    - min_repeats: minimum number of repeats
    - max_mismatches: maximum allowed mismatches per repeat unit

    Returns:
    List of dictionaries containing imperfect repeat information
    """
    try:
        sequence = str(seq).upper()
        imperfect_repeats = []

        if len(sequence) < min_period * min_repeats:
            return []

        for period in range(min_period, min(max_period + 1, len(sequence) // min_repeats + 1)):
            for start in range(len(sequence) - period * min_repeats + 1):
                pattern = sequence[start:start + period]

                # Handle patterns with ambiguous bases (N) as wildcards

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

        # Remove overlaps and sort
        imperfect_repeats = _remove_overlapping_repeats(imperfect_repeats)
        imperfect_repeats.sort(key=lambda x: (-x['total_length'], x['start']))

        return imperfect_repeats[:20]

    except Exception as e:
        logging.warning(f"Imperfect repeat analysis failed: {e}")
        return []

def find_microsatellites(seq: Seq, min_repeats=4) -> List[Dict]:
    """
    Find microsatellite repeats (1-6 bp repeat units)

    Parameters:
    - seq: BioPython Seq object
    - min_repeats: minimum number of repeats

    Returns:
    List of microsatellite repeat information
    """
    try:
        sequence = str(seq).upper()
        microsatellites = []

        # Common microsatellite patterns
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
                # Find all occurrences using regex
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

        # Remove overlaps and sort
        microsatellites = _remove_overlapping_repeats(microsatellites)
        microsatellites.sort(key=lambda x: (-x['total_length'], x['start']))

        return microsatellites[:30]

    except Exception as e:
        logging.warning(f"Microsatellite analysis failed: {e}")
        return []

def _classify_repeat_type(pattern: str, period: int) -> str:
    """Classify the type of repeat based on pattern and period"""
    if period == 1:
        return f"Homopolymer ({pattern})"
    elif period == 2:
        return f"Dinucleotide ({pattern})"
    elif period == 3:
        return f"Trinucleotide ({pattern})"
    elif period == 4:
        return f"Tetranucleotide ({pattern})"
    elif period == 5:
        return f"Pentanucleotide ({pattern})"
    elif period == 6:
        return f"Hexanucleotide ({pattern})"
    else:
        return f"Complex ({period}-mer)"

def _remove_overlapping_repeats(repeats: List[Dict]) -> List[Dict]:
    """Remove overlapping repeats, keeping the longest ones"""
    if not repeats:
        return []

    # Sort by length (descending) to prioritize longer repeats
    sorted_repeats = sorted(repeats, key=lambda x: -x['total_length'])
    non_overlapping = []

    for repeat in sorted_repeats:
        overlap = False
        for existing in non_overlapping:
            # Check for overlap
            if not (repeat['end'] <= existing['start'] or repeat['start'] >= existing['end']):
                overlap = True
                break

        if not overlap:
            non_overlapping.append(repeat)

    return non_overlapping

def analyze_repeat_distribution(seq: Seq) -> Dict:
    """
    Analyze the overall distribution of repeats in a sequence

    Parameters:
    - seq: BioPython Seq object

    Returns:
    Dictionary with repeat distribution statistics
    """
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

        # Count repeat types
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
            'repeat_density': len(all_repeats) / (seq_length / 1000),  # per kb
            'most_common_period': most_common_period,
            'repeat_types': repeat_types,
            'microsatellites_found': len(microsatellites)
        }

    except Exception as e:
        logging.warning(f"Repeat distribution analysis failed: {e}")
        return {
            'total_repeats': 0,
            'repeat_coverage': 0.0,
            'average_repeat_length': 0.0,
            'repeat_density': 0.0,
            'most_common_period': 0,
            'repeat_types': {},
            'microsatellites_found': 0
        }
