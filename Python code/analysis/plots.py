import re
import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

def sanitize_filename(name: str) -> str:
    # Καθαρίζει το όνομα αρχείου αντικαθιστώντας μη επιτρεπτούς χαρακτήρες με '_'
    return re.sub(r'[<>:"/\\|?*]', '_', name)

def format_sequence_labels(seq_ids, max_length=12, large_dataset_mode=False):
    # Μορφοποιεί τις ετικέτες ακολουθιών για να αποφεύγεται ο συνωστισμός σε γραφήματα
    # Για μεγάλα σύνολα δεδομένων (>50 ακολουθίες), χρησιμοποιεί πιο επιθετική περικοπή
    if large_dataset_mode or len(seq_ids) > 50:
        # Εμφανίζει μόνο κάθε n-οστή ετικέτα για αποφυγή συνωστισμού
        step = max(1, len(seq_ids) // 20)  # εμφάνιση το πολύ 20 ετικετών
        formatted_labels = []
        for i, seq_id in enumerate(seq_ids):
            if i % step == 0:
                if len(seq_id) <= max_length:
                    formatted_labels.append(seq_id)
                else:
                    # Περισσότερο επιθετική περικοπή
                    if max_length > 4:
                        truncated = seq_id[:max_length-2] + '..'
                    else:
                        truncated = seq_id[:max_length]
                    formatted_labels.append(truncated)
            else:
                formatted_labels.append('')  # κενή ετικέτα για διάστημα
        return formatted_labels

    # Κανονική λειτουργία για μικρότερα σύνολα δεδομένων
    formatted_labels = []
    for seq_id in seq_ids:
        if len(seq_id) <= max_length:
            formatted_labels.append(seq_id)
        else:
            #διατηρεί αρχή και τέλος όπου δυνατόν
            if max_length > 6:
                mid_point = max_length - 3
                truncated = seq_id[:mid_point] + '...'
            else:
                truncated = seq_id[:max_length]
            formatted_labels.append(truncated)

    return formatted_labels

def create_individual_sequence_plots(results, output_dir):
    """Δημιουργεί ατομικά γραφήματα για κάθε ακολουθία σε ξεχωριστό υποφάκελο"""
    plots_dir = os.path.join(output_dir, 'plots', 'individual')
    os.makedirs(plots_dir, exist_ok=True)

    individual_plot_files = {}

    for seq_id, result in results.items():
        stats = result['basic_stats']

        # Υποφάκελος για κάθε ακολουθία
        safe_seq_id = sanitize_filename(seq_id)
        seq_dir = os.path.join(plots_dir, safe_seq_id)
        os.makedirs(seq_dir, exist_ok=True)

        # Δημιουργία σχήματος με πολλαπλά υπο-γραφικά για κάθε ακολουθία
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Ολοκληρωμένη ανάλυση: {seq_id}', fontsize=16, fontweight='bold')

        # 1. Πίτα σύνθεσης νουκλεοτιδίων
        nucleotides = ['A', 'T', 'G', 'C']
        counts = [stats[nt] for nt in nucleotides]
        colors = ['red', 'blue', 'green', 'orange']

        ax1.pie(counts, labels=nucleotides, autopct='%1.1f%%', colors=colors, startangle=90)
        ax1.set_title(f'Σύνθεση Νουκλεοτιδίων\nΜήκος: {stats["length"]} bp')

        # 2. Διάγραμμα ράβδων περιεκτικότητας GC vs AT
        gc_content = stats['gc_content']
        at_content = 100 - gc_content
        ax2.bar(['GC Content', 'AT Content'], [gc_content, at_content],
                color=['lightgreen', 'lightcoral'])
        ax2.set_ylabel('Ποσοστό (%)')
        ax2.set_title('Περιεκτικότητα GC vs AT')
        ax2.set_ylim(0, 100)
        for i, v in enumerate([gc_content, at_content]):
            ax2.text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')

        # 3. Κατανομή Πουρινών vs Πυριμιδινών
        purine_count = stats['A'] + stats['G']
        pyrimidine_count = stats['T'] + stats['C']
        ax3.bar(['Πουρίνες (A+G)', 'Πυριμιδίνες (T+C)'], [purine_count, pyrimidine_count],
                color=['purple', 'cyan'])
        ax3.set_ylabel('Αριθμός')
        ax3.set_title('Κατανομή Πουρινών vs Πυριμιδινών')

        # 4. Αναλυτικό διάγραμμα ράβδων νουκλεοτιδίων
        ax4.bar(nucleotides, counts, color=colors, alpha=0.7)
        ax4.set_ylabel('Αριθμός')
        ax4.set_title('Αριθμός μεμονωμένων νουκλεοτιδίων')
        ax4.set_xlabel('Νουκλεοτίδιο')
        for i, v in enumerate(counts):
            ax4.text(i, v + max(counts)*0.01, str(v), ha='center', va='bottom')

        plt.tight_layout()

        # Αποθήκευση του αρχείου ανάλυσης png στο φάκελο ακολουθίας
        individual_filepath = os.path.join(seq_dir, f'{safe_seq_id}_analysis.png')
        plt.savefig(individual_filepath, dpi=300, bbox_inches='tight')
        plt.close()

        individual_plot_files[seq_id] = individual_filepath

        # Δημιουργία γραφήματος πολυπλοκότητας αν διατίθεται
        if 'complexity' in result and result['complexity']:
            create_complexity_plot(seq_id, result['complexity'], seq_dir)

    return individual_plot_files


def create_complexity_plot(seq_id, complexity_data, seq_dir):
    """Δημιουργεί γράφημα ανάλυσης πολυπλοκότητας για μια ακολουθία"""
    if not complexity_data.get('low_complexity_regions'):
        return

    safe_seq_id = sanitize_filename(seq_id)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    fig.suptitle(f'Ανάλυση πολυπλοκότητας ακολουθίας: {seq_id}', fontsize=14, fontweight='bold')

    # 1. Οπτικοποίηση εντροπίας Shannon
    entropy = complexity_data.get('shannon_entropy', 0)
    complexity_score = complexity_data.get('complexity_score', 0)
    metrics = ['Shannon Entropy', 'Complexity Score']
    values = [entropy, complexity_score]
    max_values = [2.0, 1.0]
    x_pos = np.arange(len(metrics))
    ax1.bar(x_pos, values, color=['darkblue', 'darkgreen'], alpha=0.7)
    for i, (val, max_val) in enumerate(zip(values, max_values)):
        ax1.axhline(y=max_val, color='red', linestyle='--', alpha=0.5)
        ax1.text(i, val + 0.05, f'{val:.3f}', ha='center', va='bottom', fontweight='bold')
    ax1.set_ylabel('Βαθμολογία')
    ax1.set_title('Μετρικές πολυπλοκότητας')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(metrics)
    ax1.set_ylim(0, max(max_values) * 1.1)

    # 2. Περιοχές χαμηλής πολυπλοκότητας
    low_complexity_regions = complexity_data.get('low_complexity_regions', [])
    if low_complexity_regions:
        region_starts = [r['start'] for r in low_complexity_regions[:10]]
        region_lengths = [r['end'] - r['start'] for r in low_complexity_regions[:10]]
        region_entropies = [r['entropy'] for r in low_complexity_regions[:10]]
        ax2.bar(range(len(region_starts)), region_lengths,
                color=plt.cm.Reds([1-e/2.0 for e in region_entropies]))
        ax2.set_ylabel('Μήκος περιοχής (bp)')
        ax2.set_xlabel('Δείκτης περιοχών χαμηλής πολυπλοκότητας')
        ax2.set_title(f'Περιοχές χαμηλής πολυπλοκότητας (Σύνολο: {len(low_complexity_regions)})')
        sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds,
                                   norm=plt.Normalize(vmin=0, vmax=2.0))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax2)
        cbar.set_label('Εντροπία Shannon')
    else:
        ax2.text(0.5, 0.5, 'Δεν βρέθηκαν περιοχές χαμηλής πολυπλοκότητας',
                 ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Περιοχές χαμηλής πολυπλοκότητας')

    plt.tight_layout()

    # Αποθήκευση του αρχείου πολυπλοκότητας σε png
    complexity_filepath = os.path.join(seq_dir, f'{safe_seq_id}_complexity.png')
    plt.savefig(complexity_filepath, dpi=300, bbox_inches='tight')
    plt.close()

def create_joint_plots(results, output_dir):
    """Δημιουργεί γραφήματα σύγκρισης για όλες τις ακολουθίες μαζί"""
    plots_dir = os.path.join(output_dir, 'plots', 'joint')
    os.makedirs(plots_dir, exist_ok=True)

    joint_plot_files = {}
    joint_plot_files['gc_comparison'] = create_enhanced_gc_comparison_plot(results, plots_dir)
    joint_plot_files['length_distribution'] = create_enhanced_length_distribution_plot(results, plots_dir)
    joint_plot_files['composition_heatmap'] = create_enhanced_composition_heatmap(results, plots_dir)
    joint_plot_files['nucleotide_boxplot'] = create_enhanced_nucleotide_boxplot(results, plots_dir)
    try:
        joint_plot_files['molecular_weight'] = create_molecular_weight_plot(results, plots_dir)
    except Exception as e:
        logging.warning(f"Αποτυχία δημιουργίας γραφήματος μοριακού βάρους: {e}")
        joint_plot_files['molecular_weight'] = None
    joint_plot_files['dashboard'] = create_comprehensive_dashboard(results, plots_dir)
    joint_plot_files['correlation_matrix'] = create_correlation_matrix(results, plots_dir)
    return joint_plot_files

def generate_all_plots(results, output_dir, generate_individual=True, generate_joint=True):
    """Δημιουργεί όλα τα γραφήματα - ατομικά σε υποφακέλους και κοινά σε έναν φάκελο"""
    try:
        plot_types = []
        if generate_individual: plot_types.append("ατομικά")
        if generate_joint: plot_types.append("κοινά")

        logging.info(f"Δημιουργία {', '.join(plot_types)} γραφημάτων για {len(results)} ακολουθίες")

        individual_plots = {}
        joint_plots = {}

        if generate_individual:
            individual_plots = create_individual_sequence_plots(results, output_dir)

        if generate_joint:
            joint_plots = create_joint_plots(results, output_dir)

        all_plots = {'individual': individual_plots, 'joint': joint_plots}
        logging.info(f"Δημιουργήθηκαν {len(individual_plots)} ατομικά γραφήματα και {len(joint_plots)} κοινά γραφήματα")
        return all_plots
    except Exception as e:
        logging.error(f"Αποτυχία δημιουργίας γραφημάτων: {e}")
        raise
