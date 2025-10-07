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
    return re.sub(r'[<>:"/\\|?*]', '_', name)

def format_sequence_labels(seq_ids, max_length=12, large_dataset_mode=False):
    # For large datasets (>50 sequences), use more aggressive optimization
    if large_dataset_mode or len(seq_ids) > 50:
        # Show only every nth label to prevent overcrowding
        step = max(1, len(seq_ids) // 20)  # Show max 20 labels
        formatted_labels = []
        for i, seq_id in enumerate(seq_ids):
            if i % step == 0:
                if len(seq_id) <= max_length:
                    formatted_labels.append(seq_id)
                else:
                    # More aggressive truncation
                    if max_length > 4:
                        truncated = seq_id[:max_length-2] + '..'
                    else:
                        truncated = seq_id[:max_length]
                    formatted_labels.append(truncated)
            else:
                formatted_labels.append('')  # Empty label for spacing
        return formatted_labels

    # Standard mode for smaller datasets
    formatted_labels = []
    for seq_id in seq_ids:
        if len(seq_id) <= max_length:
            formatted_labels.append(seq_id)
        else:
            # Intelligently truncate - keep beginning and end if possible
            if max_length > 6:
                mid_point = max_length - 3
                truncated = seq_id[:mid_point] + '...'
            else:
                truncated = seq_id[:max_length]
            formatted_labels.append(truncated)

    return formatted_labels

def create_individual_sequence_plots(results, output_dir):
    """Create individual plots for each sequence in its own subfolder"""
    plots_dir = os.path.join(output_dir, 'plots', 'individual')
    os.makedirs(plots_dir, exist_ok=True)

    individual_plot_files = {}

    for seq_id, result in results.items():
        stats = result['basic_stats']

        # 🔹 Subfolder for each sequence
        safe_seq_id = sanitize_filename(seq_id)
        seq_dir = os.path.join(plots_dir, safe_seq_id)
        os.makedirs(seq_dir, exist_ok=True)

        # Create a figure with multiple subplots for each sequence
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Comprehensive Analysis: {seq_id}', fontsize=16, fontweight='bold')

        # 1. Nucleotide composition pie chart
        nucleotides = ['A', 'T', 'G', 'C']
        counts = [stats[nt] for nt in nucleotides]
        colors = ['red', 'blue', 'green', 'orange']

        ax1.pie(counts, labels=nucleotides, autopct='%1.1f%%', colors=colors, startangle=90)
        ax1.set_title(f'Nucleotide Composition\nLength: {stats["length"]} bp')

        # 2. GC vs AT content bar chart
        gc_content = stats['gc_content']
        at_content = 100 - gc_content
        ax2.bar(['GC Content', 'AT Content'], [gc_content, at_content],
                color=['lightgreen', 'lightcoral'])
        ax2.set_ylabel('Percentage (%)')
        ax2.set_title('GC vs AT Content')
        ax2.set_ylim(0, 100)
        for i, v in enumerate([gc_content, at_content]):
            ax2.text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')

        # 3. Purine vs Pyrimidine
        purine_count = stats['A'] + stats['G']
        pyrimidine_count = stats['T'] + stats['C']
        ax3.bar(['Purines (A+G)', 'Pyrimidines (T+C)'], [purine_count, pyrimidine_count],
                color=['purple', 'cyan'])
        ax3.set_ylabel('Count')
        ax3.set_title('Purine vs Pyrimidine Distribution')

        # 4. Detailed nucleotide bar chart
        ax4.bar(nucleotides, counts, color=colors, alpha=0.7)
        ax4.set_ylabel('Count')
        ax4.set_title('Individual Nucleotide Counts')
        ax4.set_xlabel('Nucleotide')
        for i, v in enumerate(counts):
            ax4.text(i, v + max(counts)*0.01, str(v), ha='center', va='bottom')

        plt.tight_layout()

        # 🔹 Save <seqid>_analysis.png in seq folder
        individual_filepath = os.path.join(seq_dir, f'{safe_seq_id}_analysis.png')
        plt.savefig(individual_filepath, dpi=300, bbox_inches='tight')
        plt.close()

        individual_plot_files[seq_id] = individual_filepath

        # 🔹 Complexity plot (if available)
        if 'complexity' in result and result['complexity']:
            create_complexity_plot(seq_id, result['complexity'], seq_dir)

    return individual_plot_files


def create_complexity_plot(seq_id, complexity_data, seq_dir):
    """Create complexity analysis plot for a sequence"""
    if not complexity_data.get('low_complexity_regions'):
        return

    safe_seq_id = sanitize_filename(seq_id)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    fig.suptitle(f'Sequence Complexity Analysis: {seq_id}', fontsize=14, fontweight='bold')

    # Plot 1: Shannon entropy visualization
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
    ax1.set_ylabel('Score')
    ax1.set_title('Complexity Metrics')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(metrics)
    ax1.set_ylim(0, max(max_values) * 1.1)

    # Plot 2: Low complexity regions
    low_complexity_regions = complexity_data.get('low_complexity_regions', [])
    if low_complexity_regions:
        region_starts = [r['start'] for r in low_complexity_regions[:10]]
        region_lengths = [r['end'] - r['start'] for r in low_complexity_regions[:10]]
        region_entropies = [r['entropy'] for r in low_complexity_regions[:10]]
        ax2.bar(range(len(region_starts)), region_lengths,
                color=plt.cm.Reds([1-e/2.0 for e in region_entropies]))
        ax2.set_ylabel('Region Length (bp)')
        ax2.set_xlabel('Low Complexity Region Index')
        ax2.set_title(f'Low Complexity Regions (Total: {len(low_complexity_regions)})')
        sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds,
                                   norm=plt.Normalize(vmin=0, vmax=2.0))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax2)
        cbar.set_label('Shannon Entropy')
    else:
        ax2.text(0.5, 0.5, 'No low complexity regions found',
                 ha='center', va='center', transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Low Complexity Regions')

    plt.tight_layout()

    # 🔹 Save <seqid>_complexity.png in seq folder
    complexity_filepath = os.path.join(seq_dir, f'{safe_seq_id}_complexity.png')
    plt.savefig(complexity_filepath, dpi=300, bbox_inches='tight')
    plt.close()

def create_joint_plots(results, output_dir):
    """Create joint comparison plots for all sequences"""
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
        logging.warning(f"Molecular weight plot generation failed: {e}")
        joint_plot_files['molecular_weight'] = None
    joint_plot_files['dashboard'] = create_comprehensive_dashboard(results, plots_dir)
    joint_plot_files['correlation_matrix'] = create_correlation_matrix(results, plots_dir)
    return joint_plot_files


def create_enhanced_gc_comparison_plot(results, plots_dir):
    """Create enhanced GC content comparison plot"""
    seq_ids = list(results.keys())
    gc_contents = [results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]
    lengths = [results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]

    # Adaptive figure sizing for large datasets
    num_sequences = len(seq_ids)
    if num_sequences > 50:
        fig_width = max(20, num_sequences * 0.15)  # Scale width with sequence count
        fig_height = 10
    else:
        fig_width = 15
        fig_height = 6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    # Bar plot with length-based coloring
    colors = plt.cm.viridis([l/max(lengths) for l in lengths])
    bars = ax1.bar(range(len(seq_ids)), gc_contents, color=colors, alpha=0.8)

    ax1.set_xlabel('Sequences')
    ax1.set_ylabel('GC Content (%)')
    ax1.set_title('GC Content Comparison Across Sequences')
    ax1.set_xticks(range(len(seq_ids)))
    large_dataset = len(seq_ids) > 50
    formatted_labels = format_sequence_labels(seq_ids, max_length=8 if large_dataset else 10, large_dataset_mode=large_dataset)
    ax1.set_xticklabels(formatted_labels, rotation=45, ha='right', fontsize=8 if large_dataset else 10)

    # Add value labels on bars only for smaller datasets
    if not large_dataset:
        for i, (bar, gc) in enumerate(zip(bars, gc_contents)):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                     f'{gc:.1f}%', ha='center', va='bottom', fontweight='bold')

    # Add reference lines for optimal ranges
    ax1.axhline(y=40, color='green', linestyle='--', alpha=0.7, label='Optimal range')
    ax1.axhline(y=60, color='green', linestyle='--', alpha=0.7)
    ax1.axhspan(40, 60, alpha=0.1, color='green')
    ax1.legend()

    # Scatter plot: GC vs Length
    marker_size = 50 if large_dataset else 100
    ax2.scatter(lengths, gc_contents, s=marker_size, c=colors, alpha=0.7, edgecolors='black')

    # Only annotate points for smaller datasets to avoid clutter
    if not large_dataset:
        formatted_labels_scatter = format_sequence_labels(seq_ids, max_length=8)
        for i, label in enumerate(formatted_labels_scatter):
            if label:  # Only annotate non-empty labels
                ax2.annotate(label, (lengths[i], gc_contents[i]),
                             xytext=(5, 5), textcoords='offset points', fontsize=9)

    ax2.set_xlabel('Sequence Length (bp)')
    ax2.set_ylabel('GC Content (%)')
    ax2.set_title('GC Content vs Sequence Length')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'gc_comparison.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def create_enhanced_length_distribution_plot(results, plots_dir):
    """Create enhanced length distribution plot"""
    lengths = [results[seq_id]['basic_stats']['length'] for seq_id in results.keys()]
    seq_ids = list(results.keys())

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # 1. Bar plot
    colors = plt.cm.Set3(np.linspace(0, 1, len(seq_ids)))
    formatted_labels = format_sequence_labels(seq_ids, max_length=10)
    bars = ax1.bar(range(len(seq_ids)), lengths, color=colors, alpha=0.8)
    ax1.set_ylabel('Length (bp)')
    ax1.set_title('Sequence Lengths')
    ax1.set_xticks(range(len(seq_ids)))
    ax1.set_xticklabels(formatted_labels, rotation=45, ha='right')

    # Add value labels
    for bar, length in zip(bars, lengths):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(lengths)*0.01,
                 f'{length:,}', ha='center', va='bottom', fontweight='bold')

    # 2. Histogram
    ax2.hist(lengths, bins=max(3, len(lengths)//2), alpha=0.7, color='skyblue', edgecolor='black')
    ax2.set_xlabel('Length (bp)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Length Distribution Histogram')
    ax2.axvline(np.mean(lengths), color='red', linestyle='--', label=f'Mean: {np.mean(lengths):.0f}')
    ax2.legend()

    # 3. Box plot
    ax3.boxplot(lengths, patch_artist=True,
                boxprops=dict(facecolor='lightblue', alpha=0.7))
    ax3.set_ylabel('Length (bp)')
    ax3.set_title('Length Distribution Box Plot')
    ax3.set_xticklabels(['All Sequences'])

    # 4. Cumulative distribution
    sorted_lengths = sorted(lengths)
    cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
    ax4.plot(sorted_lengths, cumulative, marker='o', linewidth=2, markersize=6)
    ax4.set_xlabel('Length (bp)')
    ax4.set_ylabel('Cumulative Probability')
    ax4.set_title('Cumulative Length Distribution')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'length_distribution.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def create_enhanced_composition_heatmap(results, plots_dir):
    """Create enhanced nucleotide composition heatmap"""
    seq_ids = list(results.keys())
    nucleotides = ['A', 'T', 'G', 'C']

    # Create composition matrix (percentages)
    composition_data = []
    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        total = stats['length']
        row = [stats[nt]/total*100 for nt in nucleotides]
        composition_data.append(row)

    composition_df = pd.DataFrame(composition_data, index=seq_ids, columns=nucleotides)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # 1. Percentage heatmap
    sns.heatmap(composition_df, annot=True, fmt='.1f', cmap='RdYlBu_r',
                ax=ax1, cbar_kws={'label': 'Percentage (%)'})
    ax1.set_title('Nucleotide Composition Heatmap (%)')
    ax1.set_ylabel('Sequences')

    # 2. Count heatmap
    count_data = []
    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        row = [stats[nt] for nt in nucleotides]
        count_data.append(row)

    count_df = pd.DataFrame(count_data, index=seq_ids, columns=nucleotides)
    sns.heatmap(count_df, annot=True, fmt='d', cmap='viridis',
                ax=ax2, cbar_kws={'label': 'Count'})
    ax2.set_title('Nucleotide Count Heatmap')
    ax2.set_ylabel('Sequences')

    # 3. GC/AT ratio heatmap
    gc_at_data = []
    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        gc_content = stats['gc_content']
        at_content = 100 - gc_content
        gc_at_data.append([gc_content, at_content])

    gc_at_df = pd.DataFrame(gc_at_data, index=seq_ids, columns=['GC%', 'AT%'])
    sns.heatmap(gc_at_df, annot=True, fmt='.1f', cmap='RdBu_r',
                ax=ax3, cbar_kws={'label': 'Percentage (%)'})
    ax3.set_title('GC vs AT Content Heatmap')
    ax3.set_ylabel('Sequences')

    # 4. Purine/Pyrimidine heatmap
    purine_pyrimidine_data = []
    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        total = stats['length']
        purine_pct = (stats['A'] + stats['G']) / total * 100
        pyrimidine_pct = (stats['T'] + stats['C']) / total * 100
        purine_pyrimidine_data.append([purine_pct, pyrimidine_pct])

    pp_df = pd.DataFrame(purine_pyrimidine_data, index=seq_ids,
                         columns=['Purines %', 'Pyrimidines %'])
    sns.heatmap(pp_df, annot=True, fmt='.1f', cmap='plasma',
                ax=ax4, cbar_kws={'label': 'Percentage (%)'})
    ax4.set_title('Purine vs Pyrimidine Content')
    ax4.set_ylabel('Sequences')

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'composition_heatmap.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def create_enhanced_nucleotide_boxplot(results, plots_dir):
    """Create enhanced nucleotide distribution boxplot"""
    nucleotides = ['A', 'T', 'G', 'C']

    # Prepare data for boxplot
    nucleotide_data = {nt: [] for nt in nucleotides}
    nucleotide_percentages = {nt: [] for nt in nucleotides}

    for seq_id in results.keys():
        stats = results[seq_id]['basic_stats']
        total = stats['length']
        for nt in nucleotides:
            nucleotide_data[nt].append(stats[nt])
            nucleotide_percentages[nt].append(stats[nt]/total*100)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # 1. Count boxplot
    data_for_boxplot = [nucleotide_data[nt] for nt in nucleotides]
    box1 = ax1.boxplot(data_for_boxplot, labels=nucleotides, patch_artist=True)
    colors = ['red', 'blue', 'green', 'orange']
    for patch, color in zip(box1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax1.set_ylabel('Count')
    ax1.set_title('Nucleotide Count Distribution')
    ax1.grid(True, alpha=0.3)

    # 2. Percentage boxplot
    percentage_for_boxplot = [nucleotide_percentages[nt] for nt in nucleotides]
    box2 = ax2.boxplot(percentage_for_boxplot, labels=nucleotides, patch_artist=True)
    for patch, color in zip(box2['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax2.set_ylabel('Percentage (%)')
    ax2.set_title('Nucleotide Percentage Distribution')
    ax2.grid(True, alpha=0.3)

    # 3. Violin plot for percentages
    violin_data = [nucleotide_percentages[nt] for nt in nucleotides]
    parts = ax3.violinplot(violin_data, positions=range(1, len(nucleotides)+1), showmeans=True)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)
    ax3.set_xticks(range(1, len(nucleotides)+1))
    ax3.set_xticklabels(nucleotides)
    ax3.set_ylabel('Percentage (%)')
    ax3.set_title('Nucleotide Distribution (Violin Plot)')
    ax3.grid(True, alpha=0.3)

    # 4. Summary statistics table
    ax4.axis('tight')
    ax4.axis('off')

    # Create summary statistics
    summary_data = []
    for nt in nucleotides:
        data = nucleotide_percentages[nt]
        summary_data.append([
            nt,
            f'{np.mean(data):.1f}',
            f'{np.std(data):.1f}',
            f'{np.min(data):.1f}',
            f'{np.max(data):.1f}',
            f'{np.median(data):.1f}'
        ])

    table = ax4.table(cellText=summary_data,
                      colLabels=['Nucleotide', 'Mean %', 'Std %', 'Min %', 'Max %', 'Median %'],
                      cellLoc='center',
                      loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Color the table headers
    for i in range(len(nucleotides)):
        table[(i+1, 0)].set_facecolor(colors[i])
        table[(i+1, 0)].set_alpha(0.7)

    ax4.set_title('Summary Statistics')

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'nucleotide_boxplot.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def create_molecular_weight_plot(results, plots_dir):
    """Create molecular weight comparison plot"""
    seq_ids = list(results.keys())
    mol_weights = [results[seq_id]['basic_stats']['molecular_weight'] for seq_id in seq_ids]
    lengths = [results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]
    gc_contents = [results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # 1. Molecular weight bar chart
    colors = plt.cm.plasma([mw/max(mol_weights) for mw in mol_weights])
    formatted_labels = format_sequence_labels(seq_ids, max_length=8)
    bars = ax1.bar(range(len(seq_ids)), mol_weights, color=colors, alpha=0.8)
    ax1.set_ylabel('Molecular Weight (Da)')
    ax1.set_title('Molecular Weight Comparison')
    ax1.set_xticks(range(len(seq_ids)))
    ax1.set_xticklabels(formatted_labels, rotation=45, ha='right')

    # Add value labels
    for bar, mw in zip(bars, mol_weights):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(mol_weights)*0.01,
                 f'{mw:.0f}', ha='center', va='bottom', fontweight='bold', fontsize=8)

    # 2. MW vs Length scatter
    ax2.scatter(lengths, mol_weights, s=100, c=gc_contents, cmap='viridis',
                alpha=0.7, edgecolors='black')
    formatted_labels = format_sequence_labels(seq_ids, max_length=8)
    for i, label in enumerate(formatted_labels):
        ax2.annotate(label, (lengths[i], mol_weights[i]),
                     xytext=(5, 5), textcoords='offset points', fontsize=9)

    ax2.set_xlabel('Length (bp)')
    ax2.set_ylabel('Molecular Weight (Da)')
    ax2.set_title('Molecular Weight vs Length')
    cbar = plt.colorbar(ax2.collections[0], ax=ax2)
    cbar.set_label('GC Content (%)')
    ax2.grid(True, alpha=0.3)

    # 3. MW per nucleotide
    mw_per_nt = [mw/length for mw, length in zip(mol_weights, lengths)]
    formatted_labels = format_sequence_labels(seq_ids, max_length=8)
    ax3.bar(range(len(seq_ids)), mw_per_nt, color='lightcoral', alpha=0.8)
    ax3.set_ylabel('MW per nucleotide (Da/nt)')
    ax3.set_title('Molecular Weight per Nucleotide')
    ax3.set_xticks(range(len(seq_ids)))
    ax3.set_xticklabels(formatted_labels, rotation=45, ha='right')

    # Add value labels
    for i, mw_nt in enumerate(mw_per_nt):
        ax3.text(i, mw_nt + max(mw_per_nt)*0.01,
                 f'{mw_nt:.1f}', ha='center', va='bottom', fontweight='bold')

    # 4. MW distribution pie chart
    total_mw = sum(mol_weights)
    mw_percentages = [mw/total_mw*100 for mw in mol_weights]
    formatted_labels = format_sequence_labels(seq_ids, max_length=10)

    ax4.pie(mw_percentages, labels=formatted_labels, autopct='%1.1f%%', startangle=90)
    ax4.set_title('Molecular Weight Distribution')

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'molecular_weight.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def create_comprehensive_dashboard(results, plots_dir):
    """Create a comprehensive dashboard with key metrics"""
    seq_ids = list(results.keys())
    num_sequences = len(seq_ids)
    large_dataset = num_sequences > 50

    # Calculate various metrics
    lengths = [results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]
    gc_contents = [results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]
    mol_weights = [results[seq_id]['basic_stats']['molecular_weight'] for seq_id in seq_ids]

    # AT/GC ratios
    at_gc_ratios = []
    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        at_count = stats['A'] + stats['T']
        gc_count = stats['G'] + stats['C']
        ratio = at_count / gc_count if gc_count > 0 else 0
        at_gc_ratios.append(ratio)

    # Adaptive figure sizing for large datasets
    if large_dataset:
        fig_width = max(25, num_sequences * 0.2)
        fig_height = 20
        title_fontsize = 18
        subtitle_fontsize = 12
    else:
        fig_width = 20
        fig_height = 15
        title_fontsize = 20
        subtitle_fontsize = 14

    # Create dashboard
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(4, 4, hspace=0.4, wspace=0.4)

    # Title
    fig.suptitle(f'Comprehensive Analysis Dashboard ({num_sequences} sequences)',
                 fontsize=title_fontsize, fontweight='bold')

    # 1. Length comparison (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    bars = ax1.bar(range(len(seq_ids)), lengths, color='skyblue', alpha=0.8)
    ax1.set_title('Sequence Lengths', fontweight='bold', fontsize=subtitle_fontsize)
    ax1.set_ylabel('Length (bp)')
    ax1.set_xticks(range(len(seq_ids)))
    formatted_labels = format_sequence_labels(seq_ids, max_length=6 if large_dataset else 8, large_dataset_mode=large_dataset)
    ax1.set_xticklabels(formatted_labels, rotation=45, ha='right', fontsize=8 if large_dataset else 10)

    # 2. GC content (top middle-left)
    ax2 = fig.add_subplot(gs[0, 1])
    colors = ['red' if gc < 40 or gc > 60 else 'green' for gc in gc_contents]
    ax2.bar(range(len(seq_ids)), gc_contents, color=colors, alpha=0.8)
    ax2.set_title('GC Content', fontweight='bold', fontsize=subtitle_fontsize)
    ax2.set_ylabel('GC Content (%)')
    ax2.axhspan(40, 60, alpha=0.2, color='green', label='Optimal range')
    ax2.set_xticks(range(len(seq_ids)))
    ax2.set_xticklabels(formatted_labels, rotation=45, ha='right', fontsize=8 if large_dataset else 10)

    # 3. Molecular weight (top middle-right)
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.bar(range(len(seq_ids)), mol_weights, color='orange', alpha=0.8)
    ax3.set_title('Molecular Weight', fontweight='bold', fontsize=subtitle_fontsize)
    ax3.set_ylabel('MW (Da)')
    ax3.set_xticks(range(len(seq_ids)))
    ax3.set_xticklabels(formatted_labels, rotation=45, ha='right', fontsize=8 if large_dataset else 10)

    # 4. AT/GC ratio (top right)
    ax4 = fig.add_subplot(gs[0, 3])
    colors = ['red' if ratio < 0.5 or ratio > 2.0 else 'blue' for ratio in at_gc_ratios]
    ax4.bar(range(len(seq_ids)), at_gc_ratios, color=colors, alpha=0.8)
    ax4.set_title('AT/GC Ratio', fontweight='bold', fontsize=subtitle_fontsize)
    ax4.set_ylabel('Ratio')
    ax4.axhspan(0.5, 2.0, alpha=0.2, color='blue', label='Normal range')
    ax4.set_xticks(range(len(seq_ids)))
    ax4.set_xticklabels(formatted_labels, rotation=45, ha='right', fontsize=8 if large_dataset else 10)

    # 5. Overall nucleotide composition (middle left, spanning 2 columns)
    ax5 = fig.add_subplot(gs[1, :2])
    nucleotides = ['A', 'T', 'G', 'C']
    total_nucleotides = {nt: sum(results[seq_id]['basic_stats'][nt] for seq_id in seq_ids)
                         for nt in nucleotides}
    total_count = sum(total_nucleotides.values())
    percentages = [total_nucleotides[nt]/total_count*100 for nt in nucleotides]

    colors_nt = ['red', 'blue', 'green', 'orange']
    wedges, texts, autotexts = ax5.pie(percentages, labels=nucleotides, autopct='%1.1f%%',
                                       colors=colors_nt, startangle=90)
    ax5.set_title('Overall Nucleotide Composition', fontweight='bold')

    # 6. Length vs GC scatter (middle right, spanning 2 columns)
    ax6 = fig.add_subplot(gs[1, 2:])
    scatter = ax6.scatter(lengths, gc_contents, s=150, c=mol_weights,
                          cmap='plasma', alpha=0.7, edgecolors='black')
    for i, seq_id in enumerate(seq_ids):
        ax6.annotate(seq_id, (lengths[i], gc_contents[i]),
                     xytext=(5, 5), textcoords='offset points', fontsize=10)
    ax6.set_xlabel('Length (bp)')
    ax6.set_ylabel('GC Content (%)')
    ax6.set_title('Length vs GC Content', fontweight='bold')
    ax6.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax6)
    cbar.set_label('Molecular Weight (Da)')

    # 7. Statistics summary table (bottom, spanning all columns)
    ax7 = fig.add_subplot(gs[2:, :])
    ax7.axis('tight')
    ax7.axis('off')

    # Create comprehensive statistics table
    stats_data = []
    stats_data.append(['Metric', 'Min', 'Max', 'Mean', 'Median', 'Std Dev'])

    metrics = [
        ('Length (bp)', lengths),
        ('GC Content (%)', gc_contents),
        ('MW (Da)', mol_weights),
        ('AT/GC Ratio', at_gc_ratios)
    ]

    for metric_name, data in metrics:
        stats_data.append([
            metric_name,
            f'{np.min(data):.1f}',
            f'{np.max(data):.1f}',
            f'{np.mean(data):.1f}',
            f'{np.median(data):.1f}',
            f'{np.std(data):.1f}'
        ])

    table = ax7.table(cellText=stats_data[1:], colLabels=stats_data[0],
                      cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 2)

    # Style the table
    for i in range(len(stats_data[0])):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax7.set_title('Statistical Summary', fontweight='bold', pad=20)

    plt.savefig(os.path.join(plots_dir, 'dashboard.png'), dpi=300, bbox_inches='tight')
    plt.close()

    return os.path.join(plots_dir, 'dashboard.png')

def create_correlation_matrix(results, plots_dir):
    """Create correlation matrix of sequence metrics"""
    if len(results) < 3:
        return None  # Need at least 3 sequences for meaningful correlation

    # Prepare data for correlation analysis
    seq_ids = list(results.keys())
    data_dict = {'Sequence': seq_ids}

    # Basic metrics
    data_dict['Length'] = [results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]
    data_dict['GC_Content'] = [results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]
    data_dict['Molecular_Weight'] = [results[seq_id]['basic_stats']['molecular_weight'] for seq_id in seq_ids]

    # Nucleotide counts
    for nt in ['A', 'T', 'G', 'C']:
        data_dict[f'{nt}_Count'] = [results[seq_id]['basic_stats'][nt] for seq_id in seq_ids]
        data_dict[f'{nt}_Percent'] = [results[seq_id]['basic_stats'][nt]/results[seq_id]['basic_stats']['length']*100
                                      for seq_id in seq_ids]

    # Additional metrics
    data_dict['AT_GC_Ratio'] = []
    data_dict['Purine_Percent'] = []
    data_dict['Pyrimidine_Percent'] = []

    for seq_id in seq_ids:
        stats = results[seq_id]['basic_stats']
        at_count = stats['A'] + stats['T']
        gc_count = stats['G'] + stats['C']
        total = stats['length']

        data_dict['AT_GC_Ratio'].append(at_count / gc_count if gc_count > 0 else 0)
        data_dict['Purine_Percent'].append((stats['A'] + stats['G']) / total * 100)
        data_dict['Pyrimidine_Percent'].append((stats['T'] + stats['C']) / total * 100)

    # Create DataFrame and correlation matrix
    df = pd.DataFrame(data_dict)
    correlation_matrix = df.select_dtypes(include=[np.number]).corr()

    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # 1. Full correlation heatmap
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    sns.heatmap(correlation_matrix, mask=mask, annot=True, cmap='RdBu_r', center=0,
                square=True, ax=ax1, cbar_kws={'label': 'Correlation Coefficient'})
    ax1.set_title('Sequence Metrics Correlation Matrix', fontweight='bold')

    # 2. Clustermap for hierarchical clustering
    ax2.axis('off')

    # Create a separate figure for clustermap
    plt.figure(figsize=(12, 10))
    clustermap = sns.clustermap(correlation_matrix, annot=True, cmap='RdBu_r', center=0,
                                square=True, cbar_kws={'label': 'Correlation Coefficient'})
    clustermap.fig.suptitle('Hierarchical Clustering of Sequence Metrics',
                            fontweight='bold', y=1.02)

    # Save clustermap separately
    cluster_filename = os.path.join(plots_dir, 'correlation_clustermap.png')
    clustermap.savefig(cluster_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Add text summary to ax2
    ax2.text(0.1, 0.9, 'Correlation Analysis Summary:', fontweight='bold', fontsize=14,
             transform=ax2.transAxes)

    # Find strongest correlations
    corr_pairs = []
    for i in range(len(correlation_matrix.columns)):
        for j in range(i+1, len(correlation_matrix.columns)):
            corr_val = correlation_matrix.iloc[i, j]
            if abs(corr_val) > 0.5:  # Only show strong correlations
                corr_pairs.append((correlation_matrix.columns[i],
                                   correlation_matrix.columns[j], corr_val))

    # Sort by absolute correlation value
    corr_pairs.sort(key=lambda x: abs(x[2]), reverse=True)

    y_pos = 0.8
    ax2.text(0.1, y_pos, 'Strong Correlations (|r| > 0.5):', fontweight='bold',
             transform=ax2.transAxes)
    y_pos -= 0.08

    for i, (var1, var2, corr) in enumerate(corr_pairs[:10]):  # Show top 10
        if y_pos < 0.1:
            break
        color = 'red' if corr > 0 else 'blue'
        ax2.text(0.1, y_pos, f'{var1} ↔ {var2}: {corr:.3f}',
                 color=color, transform=ax2.transAxes)
        y_pos -= 0.06

    if not corr_pairs:
        ax2.text(0.1, 0.7, 'No strong correlations found (|r| > 0.5)',
                 transform=ax2.transAxes)

    plt.tight_layout()
    filename = os.path.join(plots_dir, 'correlation_matrix.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    return filename

def generate_all_plots(results, output_dir, generate_individual=True, generate_joint=True):
    """Generate all plots - individual (in subfolders) and joint (in one folder)"""
    try:
        plot_types = []
        if generate_individual: plot_types.append("individual")
        if generate_joint: plot_types.append("joint")

        logging.info(f"Generating {', '.join(plot_types)} plots for {len(results)} sequences")

        individual_plots = {}
        joint_plots = {}

        if generate_individual:
            individual_plots = create_individual_sequence_plots(results, output_dir)

        if generate_joint:
            joint_plots = create_joint_plots(results, output_dir)

        all_plots = {'individual': individual_plots, 'joint': joint_plots}
        logging.info(f"Generated {len(individual_plots)} individual plots and {len(joint_plots)} joint plots")
        return all_plots
    except Exception as e:
        logging.error(f"Plot generation failed: {e}")
        raise