
import os
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from .plots import generate_all_plots

def generate_pdf_report(results, output_dir):
    """Generate a comprehensive PDF report"""
    try:
        from reportlab.lib.pagesizes import letter, A4
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib import colors
        from reportlab.lib.units import inch

        filename = os.path.join(output_dir, "comprehensive_analysis_report.pdf")
        doc = SimpleDocTemplate(filename, pagesize=A4, topMargin=0.5*inch, bottomMargin=0.5*inch)
        styles = getSampleStyleSheet()
        story = []

        # Προσαρμοσμένα στυλ
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            spaceAfter=30,
            textColor=colors.darkblue,
            alignment=1 
        )

        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=16,
            spaceAfter=12,
            textColor=colors.darkgreen,
            spaceBefore=20
        )

         # Σελίδα τίτλου
        title = Paragraph("🧬 Professional DNA/RNA Analysis Report", title_style)
        story.append(title)
        story.append(Spacer(1, 30))

        # metadata αναφοράς
        metadata = f"""
        <para align="center">
        <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br/>
        <b>Number of sequences analyzed:</b> {len(results)}<br/>
        <b>Analysis Software:</b> Professional DNA/RNA Bioinformatics Suite v2.0
        </para>
        """
        story.append(Paragraph(metadata, styles['Normal']))
        story.append(Spacer(1, 30))

        # Summary
        story.append(Paragraph("📋 Executive Summary", heading_style))

        # Υπολογισμός βασικών στατιστικών
        lengths = [result['basic_stats']['length'] for result in results.values()]
        gc_contents = [result['basic_stats']['gc_content'] for result in results.values()]
        mol_weights = [result['basic_stats']['molecular_weight'] for result in results.values()]

        summary_text = f"""
        This report presents a comprehensive analysis of {len(results)} DNA/RNA sequences.
        <br/><br/>
        <b>Key Findings:</b><br/>
        • Total sequences analyzed: {len(results)}<br/>
        • Length range: {min(lengths):,} - {max(lengths):,} bp<br/>
        • Average length: {np.mean(lengths):.0f} bp<br/>
        • GC content range: {min(gc_contents):.2f}% - {max(gc_contents):.2f}%<br/>
        • Average GC content: {np.mean(gc_contents):.2f}%<br/>
        • Average molecular weight: {np.mean(mol_weights):.2f} Da
        """
        story.append(Paragraph(summary_text, styles['Normal']))
        story.append(Spacer(1, 20))

        # Πίνακας συνοπτικών στατιστικών
        story.append(Paragraph("📊 Summary Statistics", heading_style))
        summary_data = [["Sequence ID", "Length (bp)", "A", "T", "G", "C", "GC%", "Mol. Weight (Da)"]]

        for seq_id, result in results.items():
            stats = result['basic_stats']
            summary_data.append([
                seq_id,
                f"{stats['length']:,}",
                str(stats['A']),
                str(stats['T']),
                str(stats['G']),
                str(stats['C']),
                f"{stats['gc_content']:.2f}%",
                f"{stats['molecular_weight']:.2f}"
            ])

        summary_table = Table(summary_data)
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightblue),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('FONTSIZE', (0, 1), (-1, -1), 8),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.lightblue, colors.white])
        ]))
        story.append(summary_table)
        story.append(PageBreak())

       # Λεπτομερής ανάλυση ανά ακολουθία
        story.append(Paragraph("🔬 Detailed Sequence Analysis", heading_style))

        for seq_id, result in results.items():
            stats = result['basic_stats']

            # Sequence header
            seq_title = Paragraph(f"<b>Sequence: {seq_id}</b>", styles['Heading3'])
            story.append(seq_title)

            # statistics
            basic_stats_data = [
                ["Property", "Value", "Percentage"],
                ["Length", f"{stats['length']:,} bp", "100%"],
                ["Adenine (A)", str(stats['A']), f"{(stats['A']/stats['length']*100):.2f}%"],
                ["Thymine (T)", str(stats['T']), f"{(stats['T']/stats['length']*100):.2f}%"],
                ["Guanine (G)", str(stats['G']), f"{(stats['G']/stats['length']*100):.2f}%"],
                ["Cytosine (C)", str(stats['C']), f"{(stats['C']/stats['length']*100):.2f}%"],
                ["GC Content", f"{stats['gc_content']:.2f}%", ""],
                ["AT Content", f"{100-stats['gc_content']:.2f}%", ""],
                ["Molecular Weight", f"{stats['molecular_weight']:.2f} Da", ""]
            ]

            basic_table = Table(basic_stats_data, colWidths=[2*inch, 1.5*inch, 1*inch])
            basic_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            story.append(basic_table)

            # Επιπλέον αναλύσεις αν υπάρχουν
            if 'tandem_repeats' in result and result['tandem_repeats']:
                story.append(Spacer(1, 10))
                repeats_text = f"<b>Tandem Repeats Found:</b> {len(result['tandem_repeats'])}"
                story.append(Paragraph(repeats_text, styles['Normal']))

            if 'motifs' in result and result['motifs']:
                motifs_text = f"<b>Sequence Motifs Found:</b> {len(result['motifs'])}"
                story.append(Paragraph(motifs_text, styles['Normal']))

            if 'complexity' in result and result['complexity']:
                complexity = result['complexity']
                complexity_text = f"""
                <b>Complexity Analysis:</b><br/>
                • Shannon Entropy: {complexity.get('shannon_entropy', 'N/A')}<br/>
                • Complexity Score: {complexity.get('complexity_score', 'N/A')}<br/>
                • Low Complexity Regions: {len(complexity.get('low_complexity_regions', []))}
                """
                story.append(Paragraph(complexity_text, styles['Normal']))

            story.append(Spacer(1, 20))

        # visualizations
        plots_dir = os.path.join(output_dir, 'plots')
        joint_plots_dir = os.path.join(plots_dir, 'joint')
        individual_plots_dir = os.path.join(plots_dir, 'individual')

        if os.path.exists(plots_dir):
            story.append(PageBreak())
            story.append(Paragraph("📈 Data Visualizations", heading_style))

            if os.path.exists(joint_plots_dir):
                story.append(Paragraph("🔗 Comparative Analysis Plots", styles['Heading3']))

                joint_plot_files = [f for f in os.listdir(joint_plots_dir) if f.endswith('.png')]

                for plot_file in joint_plot_files:
                    plot_path = os.path.join(joint_plots_dir, plot_file)
                    if os.path.exists(plot_path):
                        try:
                            plot_title = os.path.splitext(plot_file)[0].replace('_', ' ').title()
                            story.append(Paragraph(f"<b>{plot_title}</b>", styles['Heading4']))
                            img = Image(plot_path, width=7*inch, height=5*inch)
                            story.append(img)
                            story.append(Spacer(1, 20))
                        except Exception as e:
                            logging.error(f"Error adding joint plot {plot_file}: {e}")

            # sequence plots
            if os.path.exists(individual_plots_dir):
                story.append(PageBreak())
                story.append(Paragraph("🔬 Individual Sequence Analysis Plots", styles['Heading3']))

                for seq_id in results.keys():
                    individual_seq_dir = os.path.join(individual_plots_dir, seq_id)
                    if os.path.exists(individual_seq_dir):
                        individual_plot_file = f'{seq_id}_analysis.png'
                        individual_plot_path = os.path.join(individual_seq_dir, individual_plot_file)
                        if os.path.exists(individual_plot_path):
                            try:
                                story.append(Paragraph(f"<b>Detailed Analysis: {seq_id}</b>", styles['Heading4']))
                                img = Image(individual_plot_path, width=7*inch, height=5*inch)
                                story.append(img)
                                story.append(Spacer(1, 15))
                            except Exception as e:
                                logging.error(f"Error adding individual plot for {seq_id}: {e}")

                        complexity_plot_file = f'{seq_id}_complexity.png'
                        complexity_plot_path = os.path.join(individual_seq_dir, complexity_plot_file)
                        if os.path.exists(complexity_plot_path):
                            try:
                                story.append(Paragraph(f"<b>Complexity Analysis: {seq_id}</b>", styles['Heading4']))
                                img = Image(complexity_plot_path, width=6*inch, height=4*inch)
                                story.append(img)
                                story.append(Spacer(1, 15))
                            except Exception as e:
                                logging.error(f"Error adding complexity plot for {seq_id}: {e}")

        # Quality Control Summary
        story.append(PageBreak())
        story.append(Paragraph("🛡️ Quality Control Assessment", heading_style))

        qc_issues = []
        for seq_id, result in results.items():
            stats = result['basic_stats']
            seq_issues = []

            if stats['gc_content'] < 20:
                seq_issues.append("Very low GC content (< 20%)")
            elif stats['gc_content'] > 80:
                seq_issues.append("Very high GC content (> 80%)")

            if stats['length'] < 50:
                seq_issues.append("Very short sequence (< 50 bp)")
            elif stats['length'] > 100000:
                seq_issues.append("Very long sequence (> 100kb)")

            total_standard = stats['A'] + stats['T'] + stats['G'] + stats['C']
            if total_standard < stats['length']:
                seq_issues.append(f"Contains {stats['length'] - total_standard} non-standard nucleotides")

            if seq_issues:
                qc_issues.append(f"<b>{seq_id}:</b> {'; '.join(seq_issues)}")
            else:
                qc_issues.append(f"<b>{seq_id}:</b> ✅ No issues detected")

        qc_text = "<br/>".join(qc_issues) if qc_issues else "All sequences passed quality control checks."
        story.append(Paragraph(qc_text, styles['Normal']))

        story.append(Spacer(1, 30))
        story.append(Paragraph("🔬 Analysis Methodology", heading_style))
        methodology_text = f"""
        <b>Analysis Parameters:</b><br/>
        • Sequence parsing: BioPython SeqIO<br/>
        • GC content calculation: Standard nucleotide counting<br/>
        • Molecular weight: Sum of individual nucleotide molecular weights<br/>
        • Complexity analysis: Shannon entropy calculation<br/>
        • Quality control: Industry standard thresholds<br/>
        • Visualizations: matplotlib and seaborn<br/><br/>
        
        <b>Software Version:</b> Professional DNA/RNA Bioinformatics Suite v2.0<br/>
        <b>Analysis Date:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        """
        story.append(Paragraph(methodology_text, styles['Normal']))

        # Τελική κατασκευή και αποθήκευση PDF
        doc.build(story)
        return filename

    except ImportError:
        logging.error("ReportLab not available for PDF generation")
        raise Exception("ReportLab library required for PDF generation")
    except Exception as e:
        logging.error(f"PDF generation failed: {e}")
        raise


def generate_excel_report(results, output_dir):
    """Generate a comprehensive Excel report with extensive data matching PDF"""
    try:
        filename = os.path.join(output_dir, "comprehensive_report.xlsx")

        with pd.ExcelWriter(filename, engine='xlsxwriter') as writer:
            workbook = writer.book

            # styling
            title_format = workbook.add_format({
                'bold': True,
                'font_size': 18,
                'font_color': 'navy',
                'align': 'center',
                'bg_color': '#E6F3FF',
                'border': 2
            })
            subtitle_format = workbook.add_format({
                'bold': True,
                'font_size': 14,
                'font_color': '#006400',
                'align': 'left',
                'bg_color': '#F0F8F0'
            })
            header_format = workbook.add_format({
                'bold': True,
                'bg_color': '#4472C4',
                'font_color': 'white',
                'border': 1,
                'align': 'center',
                'text_wrap': True
            })
            data_format = workbook.add_format({
                'border': 1,
                'align': 'center'
            })
            percentage_format = workbook.add_format({
                'num_format': '0.00%',
                'border': 1,
                'align': 'center'
            })
            number_format = workbook.add_format({
                'num_format': '#,##0.00',
                'border': 1,
                'align': 'center'
            })
            integer_format = workbook.add_format({
                'num_format': '#,##0',
                'border': 1,
                'align': 'center'
            })
            warning_format = workbook.add_format({
                'bg_color': '#FFCCCC',
                'border': 1,
                'align': 'center'
            })
            good_format = workbook.add_format({
                'bg_color': '#CCFFCC',
                'border': 1,
                'align': 'center'
            })

            # comprehensive statistics
            lengths = [result['basic_stats']['length'] for result in results.values()]
            gc_contents = [result['basic_stats']['gc_content'] for result in results.values()]
            mol_weights = [result['basic_stats']['molecular_weight'] for result in results.values()]
            a_counts = [result['basic_stats']['A'] for result in results.values()]
            t_counts = [result['basic_stats']['T'] for result in results.values()]
            g_counts = [result['basic_stats']['G'] for result in results.values()]
            c_counts = [result['basic_stats']['C'] for result in results.values()]

            summary_sheet = workbook.add_worksheet('Executive Summary')
            summary_sheet.merge_range('A1:E1', '🧬 Professional DNA/RNA Analysis Report', title_format)

            # metadata
            summary_sheet.write(3, 0, 'Report Generation Information:', subtitle_format)
            summary_sheet.write(4, 0, f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
            summary_sheet.write(5, 0, f'Number of sequences analyzed: {len(results)}')
            summary_sheet.write(6, 0, 'Analysis Software: Professional DNA/RNA Bioinformatics Suite v2.0')
            summary_sheet.write(7, 0, 'Platform: Python with BioPython, NumPy, Pandas, matplotlib')

            summary_sheet.write(9, 0, 'Key Statistical Findings:', subtitle_format)
            summary_sheet.write(10, 0, f'• Total sequences analyzed: {len(results)}')
            summary_sheet.write(11, 0, f'• Length range: {min(lengths):,} - {max(lengths):,} bp')
            summary_sheet.write(12, 0, f'• Average length: {np.mean(lengths):.0f} bp (Std: {np.std(lengths):.0f})')
            summary_sheet.write(13, 0, f'• Median length: {np.median(lengths):.0f} bp')
            summary_sheet.write(14, 0, f'• Total nucleotides analyzed: {sum(lengths):,} bp')
            summary_sheet.write(15, 0, f'• GC content range: {min(gc_contents):.2f}% - {max(gc_contents):.2f}%')
            summary_sheet.write(16, 0, f'• Average GC content: {np.mean(gc_contents):.2f}% (Std: {np.std(gc_contents):.2f}%)')
            summary_sheet.write(17, 0, f'• Average molecular weight: {np.mean(mol_weights):.2f} Da')
            summary_sheet.write(18, 0, f'• Total molecular weight: {sum(mol_weights):.2f} Da')

            summary_sheet.write(20, 0, 'Overall Nucleotide Composition:', subtitle_format)
            total_nucleotides = sum(a_counts) + sum(t_counts) + sum(g_counts) + sum(c_counts)
            summary_sheet.write(21, 0, f'• Total A: {sum(a_counts):,} ({sum(a_counts)/total_nucleotides*100:.2f}%)')
            summary_sheet.write(22, 0, f'• Total T: {sum(t_counts):,} ({sum(t_counts)/total_nucleotides*100:.2f}%)')
            summary_sheet.write(23, 0, f'• Total G: {sum(g_counts):,} ({sum(g_counts)/total_nucleotides*100:.2f}%)')
            summary_sheet.write(24, 0, f'• Total C: {sum(c_counts):,} ({sum(c_counts)/total_nucleotides*100:.2f}%)')
            summary_sheet.write(25, 0, f'• AT/GC Ratio: {(sum(a_counts) + sum(t_counts))/(sum(g_counts) + sum(c_counts)):.3f}')

            # Φύλλο Επισκόπησης Βασικών Στατιστικών
            overview_data = []
            for seq_id, result in results.items():
                stats = result['basic_stats']
                overview_data.append({
                    'Sequence_ID': seq_id,
                    'Length_bp': stats['length'],
                    'A_count': stats['A'],
                    'T_count': stats['T'],
                    'G_count': stats['G'],
                    'C_count': stats['C'],
                    'A_percentage': stats['A']/stats['length'],
                    'T_percentage': stats['T']/stats['length'],
                    'G_percentage': stats['G']/stats['length'],
                    'C_percentage': stats['C']/stats['length'],
                    'GC_Content': stats['gc_content']/100,
                    'AT_Content': (100-stats['gc_content'])/100,
                    'AT_GC_Ratio': (stats['A'] + stats['T'])/(stats['G'] + stats['C']) if (stats['G'] + stats['C']) > 0 else 0,
                    'Purine_Count_AG': stats['A'] + stats['G'],
                    'Pyrimidine_Count_TC': stats['T'] + stats['C'],
                    'Purine_Pyrimidine_Ratio': (stats['A'] + stats['G'])/(stats['T'] + stats['C']) if (stats['T'] + stats['C']) > 0 else 0,
                    'Molecular_Weight_Da': stats['molecular_weight'],
                    'Molecular_Weight_kDa': stats['molecular_weight']/1000,
                    'Tandem_Repeats': len(result.get('tandem_repeats', [])),
                    'Motifs_Found': len(result.get('motifs', [])),
                    'Complexity_Shannon_Entropy': result.get('complexity', {}).get('shannon_entropy', 0),
                    'Complexity_Score': result.get('complexity', {}).get('complexity_score', 0),
                    'Low_Complexity_Regions': len(result.get('complexity', {}).get('low_complexity_regions', []))
                })

            overview_df = pd.DataFrame(overview_data)
            overview_df.to_excel(writer, sheet_name='Basic Statistics Overview', index=False)

            # Επισκόπησης Βασικών Στατιστικών
            overview_sheet = writer.sheets['Basic Statistics Overview']
            for col_num, value in enumerate(overview_df.columns.values):
                overview_sheet.write(0, col_num, value, header_format)

            percentage_cols = [6, 7, 8, 9, 10, 11] 
            ratio_cols = [12, 15]  # AT/GC, Πουρινών/Πυριμιδινών
            weight_cols = [16, 17]  # Μοριακό Βάρος
            count_cols = [1, 2, 3, 4, 5, 13, 14, 18, 19, 22]  # Count columns

            for row in range(1, len(overview_df) + 1):
                for col in range(len(overview_df.columns)):
                    value = overview_df.iloc[row-1, col]
                    if col in percentage_cols:
                        overview_sheet.write(row, col, value, percentage_format)
                    elif col in weight_cols or col in [20, 21]: 
                        overview_sheet.write(row, col, value, number_format)
                    elif col in count_cols:
                        overview_sheet.write(row, col, value, integer_format)
                    elif col in ratio_cols:
                        if value < 0.5 or value > 2.0:
                            overview_sheet.write(row, col, value, warning_format)
                        else:
                            overview_sheet.write(row, col, value, good_format)
                    else:
                        overview_sheet.write(row, col, value, data_format)

            # Ανάλυση Νουκλεοτιδίων
            detailed_nucleotide_data = []
            for seq_id, result in results.items():
                stats = result['basic_stats']
                length = stats['length']

                dinucleotides = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']

                detailed_nucleotide_data.append({
                    'Sequence_ID': seq_id,
                    'Total_Length': length,
                    'Standard_Nucleotides': stats['A'] + stats['T'] + stats['G'] + stats['C'],
                    'Non_Standard_Nucleotides': length - (stats['A'] + stats['T'] + stats['G'] + stats['C']),
                    'A_Count': stats['A'],
                    'A_Frequency': stats['A']/length,
                    'T_U_Count': stats['T'],
                    'T_U_Frequency': stats['T']/length,
                    'G_Count': stats['G'],
                    'G_Frequency': stats['G']/length,
                    'C_Count': stats['C'],
                    'C_Frequency': stats['C']/length,
                    'Strong_Bonds_GC': stats['G'] + stats['C'],
                    'Weak_Bonds_AT': stats['A'] + stats['T'],
                    'Strong_Weak_Ratio': (stats['G'] + stats['C'])/(stats['A'] + stats['T']) if (stats['A'] + stats['T']) > 0 else 0,
                    'Melting_Temp_Estimate': 64.9 + 41 * (stats['G'] + stats['C'] - 16.4) / length if length > 13 else 0,
                    'Theoretical_Density': 1.7 if stats['gc_content'] > 50 else 1.6,  # Simplified density estimate
                    'Codon_Count_Estimate': length // 3,
                    'Reading_Frames_Possible': 6 if length > 100 else 3,
                    'Potential_ORFs': max(1, length // 300)  # Rough estimate
                })

            detailed_nucleotide_df = pd.DataFrame(detailed_nucleotide_data)
            detailed_nucleotide_df.to_excel(writer, sheet_name='Detailed Nucleotide Analysis', index=False)

            # Ανάλυση Νουκλεοτιδίων
            detailed_nucleotide_sheet = writer.sheets['Detailed Nucleotide Analysis']
            for col_num, value in enumerate(detailed_nucleotide_df.columns.values):
                detailed_nucleotide_sheet.write(0, col_num, value, header_format)

            # Αξιολόγηση Ποιοτικού Ελέγχου
            qc_data = []
            for seq_id, result in results.items():
                stats = result['basic_stats']
                issues = []
                warnings = []
                quality_score = 100

                if stats['length'] < 50:
                    issues.append("Very short sequence (< 50 bp)")
                    quality_score -= 30
                elif stats['length'] < 100:
                    warnings.append("Short sequence (< 100 bp)")
                    quality_score -= 10
                elif stats['length'] > 100000:
                    warnings.append("Very long sequence (> 100kb)")
                    quality_score -= 5

                if stats['gc_content'] < 20:
                    issues.append("Very low GC content (< 20%)")
                    quality_score -= 25
                elif stats['gc_content'] < 30:
                    warnings.append("Low GC content (< 30%)")
                    quality_score -= 10
                elif stats['gc_content'] > 80:
                    issues.append("Very high GC content (> 80%)")
                    quality_score -= 25
                elif stats['gc_content'] > 70:
                    warnings.append("High GC content (> 70%)")
                    quality_score -= 10

                # Σύνθεση Νουκλεοτιδίων
                total_standard = stats['A'] + stats['T'] + stats['G'] + stats['C']
                if total_standard < stats['length']:
                    non_standard_count = stats['length'] - total_standard
                    non_standard_percent = non_standard_count / stats['length'] * 100
                    if non_standard_percent > 5:
                        issues.append(f"High non-standard nucleotides ({non_standard_percent:.1f}%)")
                        quality_score -= 20
                    elif non_standard_percent > 1:
                        warnings.append(f"Contains non-standard nucleotides ({non_standard_percent:.1f}%)")
                        quality_score -= 5

                # Έλεγχοι Πολυπλοκότητας
                complexity = result.get('complexity', {})
                if complexity:
                    entropy = complexity.get('shannon_entropy', 0)
                    if entropy < 1.5:
                        warnings.append("Low sequence complexity")
                        quality_score -= 15
                    elif entropy < 1.0:
                        issues.append("Very low sequence complexity")
                        quality_score -= 25

                # Καθορισμός γενικής κατάστασης
                if quality_score >= 90:
                    status = "EXCELLENT"
                elif quality_score >= 75:
                    status = "GOOD"
                elif quality_score >= 60:
                    status = "ACCEPTABLE"
                elif quality_score >= 40:
                    status = "POOR"
                else:
                    status = "FAILED"

                qc_data.append({
                    'Sequence_ID': seq_id,
                    'Length_bp': stats['length'],
                    'GC_Content': stats['gc_content']/100,
                    'Non_Standard_Nucleotides': stats['length'] - total_standard,
                    'Non_Standard_Percentage': (stats['length'] - total_standard)/stats['length'],
                    'Shannon_Entropy': complexity.get('shannon_entropy', 0) if complexity else 0,
                    'Low_Complexity_Regions': len(complexity.get('low_complexity_regions', [])) if complexity else 0,
                    'Quality_Score': quality_score,
                    'Quality_Status': status,
                    'Issues_Count': len(issues),
                    'Warnings_Count': len(warnings),
                    'Issues': '; '.join(issues) if issues else 'None',
                    'Warnings': '; '.join(warnings) if warnings else 'None',
                    'Recommendations': _generate_qc_recommendations(stats, issues, warnings)
                })

            qc_df = pd.DataFrame(qc_data)
            qc_df.to_excel(writer, sheet_name='Quality Control Assessment', index=False)

            # Μορφοποίηση Φύλλου Ποιοτικού Ελέγχου με συνθήκη μορφοποίησης
            qc_sheet = writer.sheets['Quality Control Assessment']
            for col_num, value in enumerate(qc_df.columns.values):
                qc_sheet.write(0, col_num, value, header_format)

            excellent_format = workbook.add_format({'bg_color': '#00B050', 'font_color': 'white', 'border': 1, 'align': 'center'})
            good_format_qc = workbook.add_format({'bg_color': '#92D050', 'border': 1, 'align': 'center'})
            acceptable_format = workbook.add_format({'bg_color': '#FFFF00', 'border': 1, 'align': 'center'})
            poor_format = workbook.add_format({'bg_color': '#FFC000', 'border': 1, 'align': 'center'})
            failed_format = workbook.add_format({'bg_color': '#C00000', 'font_color': 'white', 'border': 1, 'align': 'center'})

            for row in range(1, len(qc_df) + 1):
                status = qc_df.iloc[row-1, 8]  
                for col in range(len(qc_df.columns)):
                    value = qc_df.iloc[row-1, col]
                    if col == 8:  
                        if status == 'EXCELLENT':
                            qc_sheet.write(row, col, value, excellent_format)
                        elif status == 'GOOD':
                            qc_sheet.write(row, col, value, good_format_qc)
                        elif status == 'ACCEPTABLE':
                            qc_sheet.write(row, col, value, acceptable_format)
                        elif status == 'POOR':
                            qc_sheet.write(row, col, value, poor_format)
                        else:
                            qc_sheet.write(row, col, value, failed_format)
                    elif col in [2, 4]: 
                        qc_sheet.write(row, col, value, percentage_format)
                    elif col in [5, 7]: 
                        qc_sheet.write(row, col, value, number_format)
                    elif col in [1, 3, 6, 9, 10]:
                        qc_sheet.write(row, col, value, integer_format)
                    else:
                        qc_sheet.write(row, col, value, data_format)

            # Φύλλο Συγκριτικής Ανάλυσης
            if len(results) > 1:
                comp_data = []
                seq_ids = list(results.keys())
                for i, seq1_id in enumerate(seq_ids):
                    stats1 = results[seq1_id]['basic_stats']
                    for j, seq2_id in enumerate(seq_ids):
                        if i < j:  # Αποφυγή διπλών συγκρίσεων
                            stats2 = results[seq2_id]['basic_stats']

                            # υπολογισμό λεπτομερών μετρικών
                            length_diff = abs(stats1['length'] - stats2['length'])
                            length_similarity = 1 - length_diff/max(stats1['length'], stats2['length'])

                            gc_diff = abs(stats1['gc_content'] - stats2['gc_content'])
                            gc_similarity = 1 - gc_diff/100

                            mw_diff = abs(stats1['molecular_weight'] - stats2['molecular_weight'])
                            mw_similarity = 1 - mw_diff/max(stats1['molecular_weight'], stats2['molecular_weight'])

                            # Ομοιότητα σύνθεσης νουκλεοτιδίων
                            comp1 = np.array([stats1['A'], stats1['T'], stats1['G'], stats1['C']]) / stats1['length']
                            comp2 = np.array([stats2['A'], stats2['T'], stats2['G'], stats2['C']]) / stats2['length']
                            composition_similarity = 1 - np.sum(np.abs(comp1 - comp2)) / 2

                            # Συνολική βαθμολογία ομοιότητας
                            overall_similarity = np.mean([length_similarity, gc_similarity, mw_similarity, composition_similarity])

                            comp_data.append({
                                'Sequence_1': seq1_id,
                                'Sequence_2': seq2_id,
                                'Length_1_bp': stats1['length'],
                                'Length_2_bp': stats2['length'],
                                'Length_Difference_bp': length_diff,
                                'Length_Similarity': length_similarity,
                                'GC_Content_1': stats1['gc_content']/100,
                                'GC_Content_2': stats2['gc_content']/100,
                                'GC_Difference': gc_diff/100,
                                'GC_Similarity': gc_similarity,
                                'MW_1_Da': stats1['molecular_weight'],
                                'MW_2_Da': stats2['molecular_weight'],
                                'MW_Difference_Da': mw_diff,
                                'MW_Similarity': mw_similarity,
                                'Composition_Similarity': composition_similarity,
                                'Overall_Similarity_Score': overall_similarity,
                                'Similarity_Category': _categorize_similarity(overall_similarity)
                            })

                if comp_data:
                    comp_df = pd.DataFrame(comp_data)
                    comp_df.to_excel(writer, sheet_name='Comparative Analysis', index=False)

                    # Μορφοποίηση φύλλου συγκριτικής ανάλυσης
                    comp_sheet = writer.sheets['Comparative Analysis']
                    for col_num, value in enumerate(comp_df.columns.values):
                        comp_sheet.write(0, col_num, value, header_format)

            # Φύλλο Χαρακτηριστικών Αλληλουχίας
            features_data = []
            for seq_id, result in results.items():
                stats = result['basic_stats']
                complexity = result.get('complexity', {})

                features_data.append({
                    'Sequence_ID': seq_id,
                    'Length_Category': _categorize_length(stats['length']),
                    'GC_Category': _categorize_gc_content(stats['gc_content']),
                    'Tandem_Repeats_Found': len(result.get('tandem_repeats', [])),
                    'Motifs_Found': len(result.get('motifs', [])),
                    'Has_Low_Complexity': len(complexity.get('low_complexity_regions', [])) > 0 if complexity else False,
                    'Complexity_Level': _categorize_complexity(complexity.get('shannon_entropy', 0)) if complexity else 'Unknown',
                    'Estimated_Gene_Count': max(1, stats['length'] // 1000),  # Rough estimate
                    'Potential_Coding_Percentage': min(100, max(10, stats['gc_content'] * 1.2)),  # Rough estimate
                    'Palindrome_Potential': 'High' if 45 <= stats['gc_content'] <= 55 else 'Medium' if 35 <= stats['gc_content'] <= 65 else 'Low',
                    'PCR_Suitability': _assess_pcr_suitability(stats),
                    'Cloning_Difficulty': _assess_cloning_difficulty(stats),
                    'Sequencing_Quality_Prediction': _predict_sequencing_quality(stats, complexity)
                })

            features_df = pd.DataFrame(features_data)
            features_df.to_excel(writer, sheet_name='Sequence Features', index=False)

            # Φύλλο Στατιστικής Σύνοψης
            stats_summary_data = []

            # Γενικά Στατιστικά
            stats_summary_data.extend([
                {'Statistic': 'Total Sequences', 'Value': len(results), 'Unit': 'count'},
                {'Statistic': 'Total Length', 'Value': sum(lengths), 'Unit': 'bp'},
                {'Statistic': 'Average Length', 'Value': np.mean(lengths), 'Unit': 'bp'},
                {'Statistic': 'Median Length', 'Value': np.median(lengths), 'Unit': 'bp'},
                {'Statistic': 'Length Standard Deviation', 'Value': np.std(lengths), 'Unit': 'bp'},
                {'Statistic': 'Minimum Length', 'Value': min(lengths), 'Unit': 'bp'},
                {'Statistic': 'Maximum Length', 'Value': max(lengths), 'Unit': 'bp'},
                {'Statistic': 'Length Range', 'Value': max(lengths) - min(lengths), 'Unit': 'bp'},
                {'Statistic': '', 'Value': '', 'Unit': ''},  # Separator
                {'Statistic': 'Average GC Content', 'Value': np.mean(gc_contents), 'Unit': '%'},
                {'Statistic': 'Median GC Content', 'Value': np.median(gc_contents), 'Unit': '%'},
                {'Statistic': 'GC Content Standard Deviation', 'Value': np.std(gc_contents), 'Unit': '%'},
                {'Statistic': 'Minimum GC Content', 'Value': min(gc_contents), 'Unit': '%'},
                {'Statistic': 'Maximum GC Content', 'Value': max(gc_contents), 'Unit': '%'},
                {'Statistic': 'GC Content Range', 'Value': max(gc_contents) - min(gc_contents), 'Unit': '%'},
                {'Statistic': '', 'Value': '', 'Unit': ''},  # Separator
                {'Statistic': 'Total Nucleotides', 'Value': total_nucleotides, 'Unit': 'count'},
                {'Statistic': 'Total A nucleotides', 'Value': sum(a_counts), 'Unit': 'count'},
                {'Statistic': 'Total T nucleotides', 'Value': sum(t_counts), 'Unit': 'count'},
                {'Statistic': 'Total G nucleotides', 'Value': sum(g_counts), 'Unit': 'count'},
                {'Statistic': 'Total C nucleotides', 'Value': sum(c_counts), 'Unit': 'count'},
                {'Statistic': 'Overall A percentage', 'Value': sum(a_counts)/total_nucleotides*100, 'Unit': '%'},
                {'Statistic': 'Overall T percentage', 'Value': sum(t_counts)/total_nucleotides*100, 'Unit': '%'},
                {'Statistic': 'Overall G percentage', 'Value': sum(g_counts)/total_nucleotides*100, 'Unit': '%'},
                {'Statistic': 'Overall C percentage', 'Value': sum(c_counts)/total_nucleotides*100, 'Unit': '%'},
            ])

            stats_summary_df = pd.DataFrame(stats_summary_data)
            stats_summary_df.to_excel(writer, sheet_name='Statistical Summary', index=False)

            # Φύλλο Μεθοδολογίας Ανάλυσης
            methodology_sheet = workbook.add_worksheet('Analysis Methodology')
            methodology_sheet.merge_range('A1:C1', 'Analysis Methodology & Technical Details', title_format)

            row = 3
            methodology_sheet.write(row, 0, 'Analysis Parameters:', subtitle_format)
            row += 1
            methodology_items = [
                '• Sequence parsing: BioPython SeqIO library',
                '• GC content calculation: Standard nucleotide counting method',
                '• Molecular weight: Sum of individual nucleotide molecular weights',
                '• Complexity analysis: Shannon entropy calculation',
                '• Quality control: Industry standard thresholds and best practices',
                '• Statistical analysis: NumPy and SciPy statistical functions',
                '• Data processing: Pandas DataFrame operations',
                '• Visualizations: matplotlib and seaborn libraries'
            ]

            for item in methodology_items:
                methodology_sheet.write(row, 0, item)
                row += 1

            row += 1
            methodology_sheet.write(row, 0, 'Quality Control Thresholds:', subtitle_format)
            row += 1
            qc_items = [
                '• Minimum length for analysis: 10 bp',
                '• Short sequence warning: < 100 bp',
                '• Long sequence warning: > 100,000 bp',
                '• Low GC content warning: < 30%',
                '• High GC content warning: > 70%',
                '• Very low GC content issue: < 20%',
                '• Very high GC content issue: > 80%',
                '• Non-standard nucleotide tolerance: < 5%',
                '• Low complexity threshold: Shannon entropy < 1.5'
            ]

            for item in qc_items:
                methodology_sheet.write(row, 0, item)
                row += 1

            row += 1
            methodology_sheet.write(row, 0, 'Software Information:', subtitle_format)
            row += 1
            software_items = [
                'Software Version: Professional DNA/RNA Bioinformatics Suite v2.0',
                f'Analysis Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
                'Platform: Python 3.8+ with scientific computing stack',
                'Core Libraries: BioPython, NumPy, Pandas, matplotlib, xlsxwriter',
                'Analysis Engine: Custom bioinformatics algorithms',
                'Report Generation: Automated Excel and PDF reporting system'
            ]

            for item in software_items:
                methodology_sheet.write(row, 0, item)
                row += 1

            # Φύλλο Αναφορών Διαγραμμάτων
            plot_refs_sheet = workbook.add_worksheet('Plot References')
            plot_refs_sheet.merge_range('A1:C1', 'Generated Plot Files Reference', title_format)

            plot_info_data = [
                ['Plot Type', 'Filename', 'Description'],
                ['Dashboard', 'dashboard.png', 'Comprehensive overview of all sequences'],
                ['GC Comparison', 'gc_comparison.png', 'GC content comparison across sequences'],
                ['Length Distribution', 'length_distribution.png', 'Sequence length analysis'],
                ['Composition Heatmap', 'composition_heatmap.png', 'Nucleotide composition visualization'],
                ['Nucleotide Boxplot', 'nucleotide_boxplot.png', 'Distribution analysis of nucleotides'],
                ['Molecular Weight', 'molecular_weight.png', 'Molecular weight comparison'],
                ['Correlation Matrix', 'correlation_matrix.png', 'Correlation analysis of metrics'],
                ['', '', ''],  
                ['Individual Plots', '', ''],
            ]

            # Προσθήκη αναφορών ξεχωριστών διαγραμμάτων
            for seq_id in results.keys():
                plot_info_data.append([f'{seq_id} Analysis', f'individual/{seq_id}_analysis.png', f'Detailed analysis for {seq_id}'])
                plot_info_data.append([f'{seq_id} Complexity', f'individual/{seq_id}_complexity.png', f'Complexity analysis for {seq_id}'])

            # Εγγραφή αναφορών διαγραμμάτων
            for row_num, row_data in enumerate(plot_info_data):
                for col_num, cell_data in enumerate(row_data):
                    if row_num == 0:  
                        plot_refs_sheet.write(row_num, col_num, cell_data, header_format)
                    elif row_data[0] == '' and row_data[1] == '':
                        plot_refs_sheet.write(row_num, col_num, cell_data, subtitle_format)
                    else:
                        plot_refs_sheet.write(row_num, col_num, cell_data, data_format)

            # Αυτόματη προσαρμογή στηλών
            for sheet_name in writer.sheets:
                worksheet = writer.sheets[sheet_name]
                worksheet.set_column('A:A', 25)
                worksheet.set_column('B:Z', 15)  
                worksheet.set_column('AA:ZZ', 12)

        return filename

    except Exception as e:
        logging.error(f"Excel generation failed: {e}")
        raise


def _generate_qc_recommendations(stats, issues, warnings):
    """Generate QC recommendations based on sequence statistics"""
    recommendations = []

    if stats['length'] < 100:
        recommendations.append("Consider sequence extension or validation")
    if stats['gc_content'] < 30:
        recommendations.append("Check for AT-rich regions or contamination")
    elif stats['gc_content'] > 70:
        recommendations.append("Verify GC-rich content is expected")

    total_standard = stats['A'] + stats['T'] + stats['G'] + stats['C']
    if total_standard < stats['length']:
        recommendations.append("Review non-standard nucleotides")

    if not recommendations:
        recommendations.append("Sequence meets quality standards")

    return '; '.join(recommendations)


def _categorize_similarity(similarity_score):
    """Categorize similarity scores"""
    if similarity_score >= 0.95:
        return "Very High"
    elif similarity_score >= 0.85:
        return "High"
    elif similarity_score >= 0.70:
        return "Moderate"
    elif similarity_score >= 0.50:
        return "Low"
    else:
        return "Very Low"


def _categorize_length(length):
    """Categorize sequence length"""
    if length < 100:
        return "Very Short"
    elif length < 1000:
        return "Short"
    elif length < 10000:
        return "Medium"
    elif length < 100000:
        return "Long"
    else:
        return "Very Long"


def _categorize_gc_content(gc_content):
    """Categorize GC content"""
    if gc_content < 30:
        return "AT-rich"
    elif gc_content < 45:
        return "Low GC"
    elif gc_content < 55:
        return "Balanced"
    elif gc_content < 70:
        return "GC-rich"
    else:
        return "Very GC-rich"


def _categorize_complexity(entropy):
    """Categorize sequence complexity based on Shannon entropy"""
    if entropy < 1.0:
        return "Very Low"
    elif entropy < 1.5:
        return "Low"
    elif entropy < 1.8:
        return "Medium"
    elif entropy < 2.0:
        return "High"
    else:
        return "Very High"


def _assess_pcr_suitability(stats):
    """Assess PCR suitability based on sequence characteristics"""
    gc_content = stats['gc_content']
    length = stats['length']

    if 40 <= gc_content <= 60 and length >= 100:
        return "Excellent"
    elif 30 <= gc_content <= 70 and length >= 50:
        return "Good"
    elif 20 <= gc_content <= 80:
        return "Moderate"
    else:
        return "Poor"


def _assess_cloning_difficulty(stats):
    """Assess cloning difficulty"""
    gc_content = stats['gc_content']
    length = stats['length']

    if gc_content > 75 or gc_content < 25:
        return "High"
    elif length > 10000:
        return "High"
    elif 30 <= gc_content <= 70 and length < 5000:
        return "Low"
    else:
        return "Moderate"


def _predict_sequencing_quality(stats, complexity):
    """Predict sequencing quality"""
    gc_content = stats['gc_content']
    entropy = complexity.get('shannon_entropy', 2.0) if complexity else 2.0

    if 40 <= gc_content <= 60 and entropy > 1.5:
        return "High"
    elif 30 <= gc_content <= 70 and entropy > 1.0:
        return "Good"
    else:
        return "Moderate"


def create_comprehensive_report(results, output_dir, *args, **kwargs):
    """Create comprehensive report with plots based on toggle settings - flexible signature"""
    try:
        if len(args) >= 2:
            generate_individual, generate_joint = args[:2]
        else:
            generate_individual = kwargs.get('generate_individual', True)
            generate_joint = kwargs.get('generate_joint', True)

        logging.info(f"Creating comprehensive report for {len(results)} sequences")
        logging.info(f"Generate individual plots: {generate_individual}")
        logging.info(f"Generate joint plots: {generate_joint}")
        os.makedirs(output_dir, exist_ok=True)
        plot_results = generate_all_plots(results, output_dir, generate_individual, generate_joint)
        return {
            'plots_generated': plot_results,
            'output_directory': output_dir,
            'individual_plots_enabled': generate_individual,
            'joint_plots_enabled': generate_joint
        }
    except Exception as e:
        logging.error(f"Error creating comprehensive report: {e}")

        raise

