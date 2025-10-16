import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Δημιουργεί γραφικές απεικονίσεις από τα αποτελέσματα της ανάλυσης DNA
def generate_visualizations(results, output_dir):
    plot_dir = os.path.join(output_dir, 'plots')

    # Καθαρίζει παλιές εικόνες png από τον φάκελο αποθήκευσης
    if os.path.exists(plot_dir):
        for file in os.listdir(plot_dir):
            if file.endswith('.png'):
                try:
                    os.remove(os.path.join(plot_dir, file))
                except Exception:
                    pass

    os.makedirs(plot_dir, exist_ok=True)

    try:
        # Λίστα γραφημάτων προς δημιουργία με το αντίστοιχο όνομα αρχείου και τη συνάρτηση σχεδίασης
        plots = [
            ('gc_comparison.png', create_gc_comparison_plot),
            ('length_distribution.png', create_length_distribution_plot),
            ('composition_heatmap.png', create_composition_heatmap),
            ('nucleotide_boxplot.png', create_boxplot_visualization)
        ]

        # Εκτέλεση σχεδίασης κάθε γραφήματος με κλείσιμο του plot μετά από κάθε αποθήκευση
        for filename, plot_func in plots:
            plot_func(results, plot_dir)
            plt.close('all')

    except Exception as e:
        logging.error(f"Visualization generation failed: {e}")

# Δημιουργεί ράβδους σύγκρισης GC περιεκτικότητας ανά ακολουθία
def create_gc_comparison_plot(results, plot_dir):
    seq_ids, gc_contents = zip(*[(seq_id, result['basic_stats']['gc_content']) for seq_id, result in results.items()])
    plt.figure(figsize=(10, 6))
    plt.bar(seq_ids, gc_contents)
    plt.xlabel('Sequences')
    plt.ylabel('GC Content (%)')
    plt.title('GC Content Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gc_comparison.png'), dpi=300)
    plt.close()

# Δημιουργεί ιστογράφημα κατανομής μήκους ακολουθίας
def create_length_distribution_plot(results, plot_dir):
    lengths = [result['basic_stats']['length'] for result in results.values()]
    plt.figure(figsize=(8, 6))
    plt.hist(lengths, bins=20)
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'length_distribution.png'), dpi=300)
    plt.close()

# Δημιουργεί θερμικό χάρτη με τις περιεκτικότητες των νουκλεοτιδίων
def create_composition_heatmap(results, plot_dir):
    data = [[result['basic_stats'][base] for base in 'ATGC'] for result in results.values()]
    df = pd.DataFrame(data, columns=['A', 'T', 'G', 'C'], index=list(results.keys()))
    plt.figure(figsize=(8, 6))
    sns.heatmap(df, annot=True, fmt='d', cmap='YlOrRd')
    plt.title('Nucleotide Composition Heatmap')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'composition_heatmap.png'), dpi=300)
    plt.close()

# Δημιουργεί boxplot για τη διασπορά ποσοστών των νουκλεοτιδίων στις ακολουθίες
def create_boxplot_visualization(results, plot_dir):
    plot_data = []
    for seq_id, result in results.items():
        stats = result['basic_stats']
        total_length = stats['length']
        if total_length > 0:
            for nucleotide in 'ATGC':
                plot_data.append({
                    'Nucleotide': nucleotide,
                    'Percentage': stats[nucleotide] / total_length * 100
                })

    df_boxplot = pd.DataFrame(plot_data)
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_boxplot, x='Nucleotide', y='Percentage')
    plt.title('Nucleotide Composition Distribution (Boxplot)')
    plt.xlabel('Nucleotide Type')
    plt.ylabel('Percentage (%)')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'nucleotide_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
