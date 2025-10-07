import os
import time
import re
import threading
import logging
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import Image, ImageTk
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import ttkbootstrap as tb
from Bio import SeqIO
from Bio.Seq import Seq
from analysis.sequences import comprehensive_analysis
from analysis.reports import create_comprehensive_report, generate_pdf_report, generate_excel_report
from utils.memory import low_memory_mode
from utils.files import save_session, load_session
from gui.dialogs import sequence_converter_window, ask_string, ask_integer

class EnhancedDNAToolGUI:
    def __init__(self, root):
        self.root = root
        root.title('Professional RNA Analysis Bioinformatics')
        root.geometry('1400x900')
        self.session = {}
        self.results = {}
        self.setup_widgets()

    def setup_widgets(self):
        nb = tb.Notebook(self.root)
        nb.pack(fill='both', expand=True)

        self.frames = {}
        tabs = [('main', 'Main Analysis'), ('advanced', 'Advanced Analysis'),
                ('qc', 'Quality Control'), ('compare', 'Sequence Comparison'), ('results', 'Results')]

        for key, text in tabs:
            self.frames[key] = tb.Frame(nb)
            nb.add(self.frames[key], text=text)

        self.setup_main_tab()
        self.setup_advanced_tab()
        self.setup_qc_tab()
        self.setup_comparison_tab()
        self.setup_results_tab()
        self.setup_menu()

    def setup_main_tab(self):
        tb.Label(self.frames['main'], text='Select FASTA/GenBank files:', font=('Arial', 12, 'bold')).pack(pady=6)
        self.file_entry = tb.Entry(self.frames['main'], width=120)
        self.file_entry.pack(pady=3)

        btn_frame = tb.Frame(self.frames['main'])
        btn_frame.pack(pady=5)
        tb.Button(btn_frame, text='Browse Files', bootstyle='primary', command=self.browse_files).pack(side='left', padx=5)

        tb.Label(self.frames['main'], text='Output folder:', font=('Arial', 12, 'bold')).pack(pady=(10, 6))
        self.output_entry = tb.Entry(self.frames['main'], width=120)
        self.output_entry.pack(pady=3)
        tb.Button(self.frames['main'], text='Select Output Folder', bootstyle='primary', command=self.select_output_folder).pack(pady=3)

        self.progress = tb.Progressbar(self.frames['main'], length=600, mode='determinate')
        self.progress.pack(pady=10)
        self.status_label = tb.Label(self.frames['main'], text='Ready', font=('Arial', 10))
        self.status_label.pack(pady=5)

        toggles_frame = tb.Frame(self.frames['main'])
        toggles_frame.pack(pady=10)

        tb.Label(toggles_frame, text='Plot Generation Options:', font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)

        self.individual_plots_var = tk.BooleanVar(value=True)
        individual_check = tb.Checkbutton(toggles_frame, text='Generate Individual Plots', variable=self.individual_plots_var, bootstyle='success')
        individual_check.grid(row=1, column=0, padx=20, pady=5, sticky='w')

        self.joint_plots_var = tk.BooleanVar(value=True)
        joint_check = tb.Checkbutton(toggles_frame, text='Generate Joint Plots', variable=self.joint_plots_var, bootstyle='success')
        joint_check.grid(row=1, column=1, padx=20, pady=5, sticky='w')

        tb.Button(self.frames['main'], text='Start Comprehensive Analysis', bootstyle='success-outline', command=self.start_analysis).pack(pady=15)

    def setup_advanced_tab(self):
        align_frame = tb.LabelFrame(self.frames['advanced'], text='Sequence Alignment', padding=10)
        align_frame.pack(pady=10, padx=20, fill='x')

        tb.Label(align_frame, text='Select two sequences to align:').pack(anchor='w')
        self.seq1_combo = ttk.Combobox(align_frame, width=50, state='readonly')
        self.seq1_combo.pack(pady=5, fill='x')
        self.seq2_combo = ttk.Combobox(align_frame, width=50, state='readonly')
        self.seq2_combo.pack(pady=5, fill='x')
        tb.Button(align_frame, text='Perform Alignment', bootstyle='info', command=self.perform_alignment).pack(pady=5)

        manual_frame = tb.LabelFrame(self.frames['advanced'], text='Manual Sequence Input', padding=10)
        manual_frame.pack(pady=10, padx=20, fill='both', expand=True)

        tb.Label(manual_frame, text='Paste sequence:').pack(anchor='w')
        self.manual_seq_text = tk.Text(manual_frame, height=8, width=80)
        self.manual_seq_text.pack(pady=5, fill='both', expand=True)

        manual_btn_frame = tb.Frame(manual_frame)
        manual_btn_frame.pack(pady=5)

        tb.Button(manual_btn_frame, text='Analyze Single Sequence', bootstyle='warning', command=self.analyze_manual_sequence).pack(side='left', padx=5)
        tb.Button(manual_btn_frame, text='Clear', bootstyle='secondary', command=lambda: self.manual_seq_text.delete(1.0, tk.END)).pack(side='left', padx=5)

    def setup_qc_tab(self):
        tb.Label(self.frames['qc'], text='Sequence Quality Control', font=('Arial', 14, 'bold')).pack(pady=10)

        qc_buttons = [
            ('Check Sequence Quality', 'primary', self.check_sequence_quality),
            ('Validate Sequence Format', 'warning', self.validate_sequence_format),
            ('Calculate AT/GC Ratio', 'info', self.calculate_at_gc_ratios),
            ('Sliding Window Analysis', 'success', self.sliding_window_analysis),
            ('Analyze Sequence Complexity', 'info', self.analyze_complexity),
            ('Find Low Complexity Regions', 'secondary', self.find_low_complexity)
        ]

        for text, style, command in qc_buttons:
            tb.Button(self.frames['qc'], text=text, bootstyle=style, command=command).pack(pady=2)

        self.qc_text = tk.Text(self.frames['qc'], height=15, width=100)
        self.qc_text.pack(pady=10, padx=20, fill='both', expand=True)

    def setup_comparison_tab(self):
        tb.Label(self.frames['compare'], text='Sequence Comparison Tools', font=('Arial', 14, 'bold')).pack(pady=10)

        msa_frame = tb.LabelFrame(self.frames['compare'], text='Multiple Sequence Alignment', padding=10)
        msa_frame.pack(pady=10, padx=20, fill='x')
        tb.Button(msa_frame, text='Compare All Sequences', bootstyle='info', command=self.compare_all_sequences).pack(pady=5)

        self.comparison_text = tk.Text(self.frames['compare'], height=20, width=100)
        self.comparison_text.pack(pady=10, padx=20, fill='both', expand=True)

    def setup_results_tab(self):
        results_nb = tb.Notebook(self.frames['results'])
        results_nb.pack(fill='both', expand=True, padx=10, pady=10)

        self.summary_frame = tb.Frame(results_nb)
        results_nb.add(self.summary_frame, text='Summary')

        summary_button_frame = tb.Frame(self.summary_frame)
        summary_button_frame.pack(pady=5)
        tb.Button(summary_button_frame, text='Export Summary to CSV', bootstyle='info', command=self.export_summary_csv).pack(side='left', padx=5)

        self.summary_text = tk.Text(self.summary_frame, height=25, width=100)
        self.summary_text.pack(fill='both', expand=True, padx=5, pady=5)

        self.detailed_frame = tb.Frame(results_nb)
        results_nb.add(self.detailed_frame, text='Detailed Results')

        detailed_button_frame = tb.Frame(self.detailed_frame)
        detailed_button_frame.pack(pady=5)
        tb.Button(detailed_button_frame, text='Export Detailed Results to CSV', bootstyle='info', command=self.export_detailed_csv).pack(side='left', padx=5)

        self.detailed_text = tk.Text(self.detailed_frame, height=25, width=100)
        self.detailed_text.pack(fill='both', expand=True, padx=5, pady=5)

        self.plots_frame = tb.Frame(results_nb)
        results_nb.add(self.plots_frame, text='Visualizations')

        # Add visualization controls
        self.setup_visualization_controls()

        self.plot_images = []
        self.plot_labels = []

    def setup_menu(self):
        menubar = tk.Menu(self.root)

        file_menu = tk.Menu(menubar, tearoff=0)
        for item in [('New Analysis', self.new_analysis), ('Save Session', self.save_session),
                     ('Load Session', self.load_session), ('---', None), ('Exit', self.root.quit)]:
            if item[1]:
                file_menu.add_command(label=item[0], command=item[1])
            else:
                file_menu.add_separator()
        menubar.add_cascade(label='File', menu=file_menu)

        tools_menu = tk.Menu(menubar, tearoff=0)
        for item in [('Sequence Converter', sequence_converter_window), ('Reverse Complement', self.reverse_complement_tool),
                     ('Translation Tool', self.translation_tool), ('Generate Report', self.generate_report),
                     ('---', None), ('Batch File Processor', self.batch_file_processor),
                     ('Sequence Length Filter', self.sequence_length_filter)]:
            if item[1]:
                tools_menu.add_command(label=item[0], command=item[1])
            else:
                tools_menu.add_separator()
        menubar.add_cascade(label='Tools', menu=tools_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label='User Guide', command=self.show_help)
        help_menu.add_command(label='About', command=self.show_about)
        menubar.add_cascade(label='Help', menu=help_menu)

        self.root.config(menu=menubar)

    def browse_files(self):
        files = filedialog.askopenfilenames(
            title='Select FASTA/GenBank files',
            filetypes=[('FASTA files', '*.fasta *.fa *.fas'), ('GenBank files', '*.gb *.gbk'), ('All supported', '*.fasta *.fa *.fas *.gb *.gbk')]
        )
        if files:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, ';'.join(files))

    def select_output_folder(self):
        folder = filedialog.askdirectory(title='Select output folder')
        if folder:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, folder)

    def start_analysis(self):
        files = self.file_entry.get().split(';') if self.file_entry.get() else []
        output_dir = self.output_entry.get()

        if not files or not output_dir:
            messagebox.showwarning('Warning', 'Please select files and output directory')
            return

        thread = threading.Thread(target=self.run_analysis_thread, args=(files, output_dir))
        thread.daemon = True
        thread.start()

    def run_analysis_thread(self, files, output_dir):
        try:
            import gc
            start_time = time.time()   # ⏱ Start timer
            self.update_status("Starting analysis...")

            is_low_memory = low_memory_mode()
            if is_low_memory:
                self.update_status("Low memory mode enabled - some features disabled")

            self.update_status("Counting sequences...")
            total_sequences = sum(
                1 for file_path in files
                for _ in SeqIO.parse(
                    file_path,
                    'genbank' if file_path.lower().endswith(('.gb', '.gbk')) else 'fasta'
                )
            )

            if total_sequences == 0:
                messagebox.showerror("Error", "No valid sequences found in the selected files")
                return

            self.update_status(f"Found {total_sequences} sequences to analyze...")
            self.update_progress(0)

            all_results = {}
            seq_names = []
            total_processed = 0

            for seq_id, seq in self.parse_files(files):
                total_processed += 1
                self.update_status(f"Analyzing {seq_id} ({total_processed}/{total_sequences})...")

                if is_low_memory and len(seq) > 50000:
                    logging.warning(f"Skipping large sequence {seq_id} in low memory mode")
                    self.update_progress((total_processed / total_sequences) * 80)
                    continue

                try:
                    result = comprehensive_analysis(seq, seq_id)
                    all_results[seq_id] = result
                    seq_names.append(seq_id)

                    del seq, result
                    gc.collect()

                    self.update_progress((total_processed / total_sequences) * 80)

                except Exception as e:
                    logging.error(f"Failed to analyze {seq_id}: {e}")
                    self.update_progress((total_processed / total_sequences) * 80)
                    continue

            if not all_results:
                messagebox.showerror("Error", "No valid sequences found or all sequences too large")
                return

            self.root.after(0, lambda: self.update_sequence_combos(seq_names))
            self.results = all_results

            self.update_status("Creating comprehensive report and visualizations...")
            self.update_progress(85)
            try:
                generate_individual = self.individual_plots_var.get()
                generate_joint = self.joint_plots_var.get()
                report_results = create_comprehensive_report(
                    all_results, output_dir, generate_individual, generate_joint
                )
                logging.info(f"Report generation completed: {report_results}")
            except Exception as e:
                logging.error(f"Report generation failed: {e}")

            gc.collect()

            self.update_status("Finalizing results...")
            self.update_progress(95)
            self.root.after(0, lambda: self.display_results(all_results))
            self.root.after(0, lambda: self.load_plots())
            self.update_progress(100)

            # ⏱ Stop timer
            end_time = time.time()
            elapsed = end_time - start_time
            elapsed_str = f"{elapsed:.2f} seconds ({elapsed/60:.2f} minutes)"

            # Permanent logging
            logging.info(f"Analysis completed in {elapsed_str} | Sequences processed: {total_processed}/{total_sequences}")

            try:
                if output_dir:
                    run_log = os.path.join(output_dir, "analysis.log")
                    with open(run_log, "a", encoding="utf-8") as f:
                        f.write(
                            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] "
                            f"Analysis completed in {elapsed_str} | "
                            f"Sequences processed: {total_processed}/{total_sequences}\n"
                        )
            except Exception as e:
                logging.warning(f"Could not write analysis time log: {e}")

            # Update GUI
            self.update_status(f"Analysis completed in {elapsed_str}")
            messagebox.showinfo(
                "Success",
                f"Analysis completed in {elapsed_str}!\n"
                f"Sequences processed: {total_processed}/{total_sequences}\n"
                f"Results saved to {output_dir}"
            )

        except Exception as e:
            logging.exception("Analysis failed")
            self.update_progress(0)
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
            self.update_status("Analysis failed")

    def parse_files(self, files):
        for file_path in files:
            try:
                format_type = 'genbank' if file_path.lower().endswith(('.gb', '.gbk')) else 'fasta'
                for record in SeqIO.parse(file_path, format_type):
                    yield (record.id, record.seq)
            except Exception as e:
                logging.error(f"Failed to parse {file_path}: {e}")

    def update_status(self, message):
        def update():
            self.status_label.config(text=message)
            logging.info(message)
        self.root.after(0, update)

    def update_progress(self, value):
        def update():
            self.progress['value'] = value
        self.root.after(0, update)

    def update_sequence_combos(self, seq_names):
        self.seq1_combo['values'] = seq_names
        self.seq2_combo['values'] = seq_names

    def display_results(self, results):
        self.summary_text.delete(1.0, tk.END)
        self.detailed_text.delete(1.0, tk.END)

        summary = "=== ANALYSIS SUMMARY ===\n\n"
        for seq_id, result in results.items():
            stats = result['basic_stats']
            summary += f"Sequence: {seq_id}\n  Length: {stats['length']} bp\n  GC Content: {stats['gc_content']:.2f}%\n\n"
        self.summary_text.insert(tk.END, summary)

        detailed = "=== DETAILED RESULTS ===\n\n"
        for seq_id, result in results.items():
            detailed += f"{'='*50}\nSEQUENCE: {seq_id}\n{'='*50}\n\n"
            stats = result['basic_stats']
            detailed += "BASIC STATISTICS:\n"
            for key, value in stats.items():
                detailed += f"  {key}: {value}\n"
            detailed += "\n"
        self.detailed_text.insert(tk.END, detailed)

    def setup_visualization_controls(self):
        """Setup dropdown and controls for visualization mode"""
        controls_frame = tb.Frame(self.plots_frame)
        controls_frame.pack(fill='x', padx=10, pady=5)

        # Visualization mode selection
        tb.Label(controls_frame, text='Visualization Mode:', font=('Arial', 12, 'bold')).pack(side='left', padx=(0, 10))

        self.viz_mode_var = tk.StringVar(value='Joint Plots')
        self.viz_mode_combo = tb.Combobox(controls_frame, textvariable=self.viz_mode_var,
                                          values=['Joint Plots', 'Individual Plots'],
                                          state='readonly', width=20)
        self.viz_mode_combo.pack(side='left', padx=(0, 10))
        self.viz_mode_combo.bind('<<ComboboxSelected>>', self.on_viz_mode_change)

        # Individual plots folder selection (initially hidden)
        self.folder_selection_frame = tb.Frame(controls_frame)
        self.folder_selection_frame.pack(side='left', padx=(20, 0))

        self.folder_label = tb.Label(self.folder_selection_frame, text='Select Sequence:', font=('Arial', 10, 'bold'))
        self.folder_label.pack(side='left', padx=(0, 5))

        self.folder_var = tk.StringVar()
        self.folder_combo = tb.Combobox(self.folder_selection_frame, textvariable=self.folder_var,
                                        state='readonly', width=25)
        self.folder_combo.pack(side='left', padx=(0, 10))
        self.folder_combo.bind('<<ComboboxSelected>>', self.on_folder_change)

        # Navigation buttons for individual plots
        self.nav_frame = tb.Frame(self.folder_selection_frame)
        self.nav_frame.pack(side='left', padx=(10, 0))

        self.prev_button = tb.Button(self.nav_frame, text='← Previous', bootstyle='secondary',
                                     command=self.previous_folder, width=10)
        self.prev_button.pack(side='left', padx=(0, 5))

        self.next_button = tb.Button(self.nav_frame, text='Next →', bootstyle='secondary',
                                     command=self.next_folder, width=10)
        self.next_button.pack(side='left')

        # Initially hide individual plot controls
        self.folder_selection_frame.pack_forget()

        # Scrollable frame for plots
        self.plots_canvas = tb.Canvas(self.plots_frame)
        self.plots_scrollbar = tb.Scrollbar(self.plots_frame, orient='vertical', command=self.plots_canvas.yview)
        self.scrollable_plots_frame = tb.Frame(self.plots_canvas)

        self.scrollable_plots_frame.bind('<Configure>', lambda e: self.plots_canvas.configure(scrollregion=self.plots_canvas.bbox('all')))

        self.plots_canvas.create_window((0, 0), window=self.scrollable_plots_frame, anchor='nw')
        self.plots_canvas.configure(yscrollcommand=self.plots_scrollbar.set)

        self.plots_canvas.pack(side='left', fill='both', expand=True, padx=(10, 0), pady=10)
        self.plots_scrollbar.pack(side='right', fill='y', pady=10)

    def on_viz_mode_change(self, event=None):
        """Handle visualization mode change"""
        mode = self.viz_mode_var.get()
        if mode == 'Individual Plots':
            self.folder_selection_frame.pack(side='left', padx=(20, 0))
            self.update_folder_list()
        else:
            self.folder_selection_frame.pack_forget()
        self.load_plots()

    def update_folder_list(self):
        """Update the list of available individual plot folders"""
        output_dir = self.output_entry.get()
        if not output_dir:
            return

        individual_plots_dir = os.path.join(output_dir, 'plots', 'individual')
        if os.path.exists(individual_plots_dir):
            folders = [f for f in os.listdir(individual_plots_dir)
                       if os.path.isdir(os.path.join(individual_plots_dir, f))]
            folders.sort()
            self.folder_combo['values'] = folders
            if folders and not self.folder_var.get():
                self.folder_var.set(folders[0])

    def on_folder_change(self, event=None):
        """Handle individual plot folder selection change"""
        self.load_plots()

    def previous_folder(self):
        """Navigate to previous folder"""
        current_folders = list(self.folder_combo['values'])
        current_folder = self.folder_var.get()
        if current_folder in current_folders:
            current_index = current_folders.index(current_folder)
            if current_index > 0:
                self.folder_var.set(current_folders[current_index - 1])
                self.load_plots()

    def next_folder(self):
        """Navigate to next folder"""
        current_folders = list(self.folder_combo['values'])
        current_folder = self.folder_var.get()
        if current_folder in current_folders:
            current_index = current_folders.index(current_folder)
            if current_index < len(current_folders) - 1:
                self.folder_var.set(current_folders[current_index + 1])
                self.load_plots()

    def load_plots(self):
        """Load plots based on selected visualization mode"""
        # Clear existing plots
        for widget in self.scrollable_plots_frame.winfo_children():
            widget.destroy()
        self.plot_images = []
        self.plot_labels = []

        output_dir = self.output_entry.get()
        if not output_dir:
            return

        viz_mode = self.viz_mode_var.get()

        if viz_mode == 'Joint Plots':
            self.load_joint_plots(output_dir)
        else:  # Individual Plots
            self.load_individual_plots(output_dir)

    def load_joint_plots(self, output_dir):
        """Load and display joint plots"""
        joint_plots_dir = os.path.join(output_dir, 'plots', 'joint')
        if not os.path.exists(joint_plots_dir):
            # Show message if no joint plots available
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text='No joint plots available. Make sure "Generate Joint Plots" is enabled and run analysis.',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        row, col = 0, 0
        max_cols = 2  # Show 2 plots per row for better clarity

        # Priority plots for better organization
        priority_plots = ['dashboard.png', 'gc_comparison.png', 'length_distribution.png',
                          'composition_heatmap.png', 'nucleotide_boxplot.png', 'molecular_weight.png',
                          'correlation_matrix.png']

        # Load priority plots first
        for plot_file in priority_plots:
            image_path = os.path.join(joint_plots_dir, plot_file)
            if os.path.exists(image_path):
                self.load_and_display_plot(image_path, plot_file, row, col, 600, 450)  # Larger size for joint plots
                col += 1
                if col >= max_cols:
                    col = 0
                    row += 1

        # Load any remaining plots
        for plot_file in sorted(os.listdir(joint_plots_dir)):
            if plot_file.endswith('.png') and plot_file not in priority_plots:
                image_path = os.path.join(joint_plots_dir, plot_file)
                self.load_and_display_plot(image_path, plot_file, row, col, 600, 450)
                col += 1
                if col >= max_cols:
                    col = 0
                    row += 1

    def load_individual_plots(self, output_dir):
        """Load and display individual plots for selected sequence"""
        individual_plots_dir = os.path.join(output_dir, 'plots', 'individual')
        if not os.path.exists(individual_plots_dir):
            # Show message if no individual plots available
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text='No individual plots available. Make sure "Generate Individual Plots" is enabled and run analysis.',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        selected_folder = self.folder_var.get()
        if not selected_folder:
            # Show message to select a sequence
            select_msg_label = tb.Label(self.scrollable_plots_frame,
                                        text='Please select a sequence from the dropdown above to view its individual plots.',
                                        font=('Arial', 12), foreground='blue')
            select_msg_label.pack(pady=20)
            self.plot_labels.append(select_msg_label)
            return

        seq_folder_path = os.path.join(individual_plots_dir, selected_folder)
        if not os.path.exists(seq_folder_path):
            # Show error message
            error_label = tb.Label(self.scrollable_plots_frame,
                                   text=f'Plots for sequence "{selected_folder}" not found.',
                                   font=('Arial', 12), foreground='red')
            error_label.pack(pady=20)
            self.plot_labels.append(error_label)
            return

        # Show sequence info header
        header_label = tb.Label(self.scrollable_plots_frame,
                                text=f'Individual Analysis Plots for: {selected_folder}',
                                font=('Arial', 14, 'bold'), foreground='navy')
        header_label.pack(pady=(10, 20))
        self.plot_labels.append(header_label)

        # Load plots for the selected sequence
        plot_files = [f for f in os.listdir(seq_folder_path) if f.endswith('.png')]
        if not plot_files:
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text=f'No plot files found for sequence "{selected_folder}".',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        # Display plots in a single column for individual sequence (easier to view)
        for i, plot_file in enumerate(sorted(plot_files)):
            image_path = os.path.join(seq_folder_path, plot_file)
            self.load_and_display_plot(image_path, plot_file, i, 0, 800, 600, single_column=True)

    def load_and_display_plot(self, image_path, plot_name, row, col, width=400, height=300, single_column=False):
        """Load and display a single plot image"""
        try:
            img = Image.open(image_path).resize((width, height))
            photo = ImageTk.PhotoImage(img)
            self.plot_images.append(photo)

            # Create frame for plot with title
            plot_frame = tb.Frame(self.scrollable_plots_frame)
            if single_column:
                plot_frame.pack(pady=10, padx=20, fill='x')
            else:
                plot_frame.grid(row=row, column=col, padx=10, pady=10, sticky='nsew')

            # Add plot title
            title = plot_name.replace('.png', '').replace('_', ' ').title()
            title_label = tb.Label(plot_frame, text=title, font=('Arial', 11, 'bold'))
            title_label.pack(pady=(0, 5))

            # Add plot image
            img_label = tb.Label(plot_frame, image=photo)
            img_label.pack()

            self.plot_labels.extend([plot_frame, title_label, img_label])

        except Exception as e:
            logging.error(f"Error loading plot {plot_name}: {e}")

    def perform_alignment(self):
        seq1_name = self.seq1_combo.get()
        seq2_name = self.seq2_combo.get()

        if not seq1_name or not seq2_name or seq1_name == seq2_name:
            messagebox.showwarning("Warning", "Please select two different sequences")
            return

        from analysis.alignments import sequence_alignment

        # Ask user for method
        method = ask_string("Alignment Method", "Choose method (NW = Needleman-Wunsch, SW = Smith-Waterman):", initial="NW")
        if not method:
            return
        method = method.strip().lower()
        if method not in ("nw", "sw"):
            messagebox.showerror("Error", "Invalid method. Please enter NW or SW.")
            return

        # Retrieve sequences from original files
        files = self.file_entry.get().split(";")
        seqs = {}
        for seq_id, seq in self.parse_files(files):
            if seq_id in (seq1_name, seq2_name):
                seqs[seq_id] = seq
            if len(seqs) == 2:
                break

        if seq1_name not in seqs or seq2_name not in seqs:
            messagebox.showerror("Error", "Could not retrieve sequences for alignment")
            return

        seq1 = seqs[seq1_name]
        seq2 = seqs[seq2_name]

        # Run alignment
        alignment_result = sequence_alignment(seq1, seq2, method=method)

        if not alignment_result:
            messagebox.showerror("Error", "Sequence alignment failed.")
            return

        # Show results
        window = tk.Toplevel(self.root)
        window.title(f"Sequence Alignment Result ({method.upper()})")
        window.geometry("800x600")

        text_widget = tk.Text(window, wrap=tk.WORD)
        text_widget.pack(fill="both", expand=True, padx=10, pady=10)

        content = f"Method: {'Needleman–Wunsch' if method=='nw' else 'Smith–Waterman'}\n"
        content += f"Alignment Score: {alignment_result['score']}\n"
        content += f"Identity: {alignment_result['identity']:.2f}%\n"
        content += f"Gaps: {alignment_result['gaps']}\n\n"
        content += "Alignment:\n" + alignment_result['alignment'] + "\n"

        text_widget.insert(tk.END, content)

    def analyze_manual_sequence(self):
        sequence_text = self.manual_seq_text.get(1.0, tk.END).strip()
        if not sequence_text:
            messagebox.showwarning("Warning", "Please enter a sequence")
            return

        sequence = re.sub(r'[^ATGCN]', '', sequence_text.upper())
        if not sequence:
            messagebox.showwarning("Warning", "No valid DNA/RNA sequence found")
            return

        seq_obj = Seq(sequence)
        result = comprehensive_analysis(seq_obj, "Manual_Input")

        window = tk.Toplevel(self.root)
        window.title("Manual Sequence Analysis")
        window.geometry("600x500")


        text_widget = tk.Text(window, wrap=tk.WORD)
        text_widget.pack(fill="both", expand=True, padx=10, pady=10)

        content = f"""Manual Sequence Analysis Results
=================================

Length: {result['basic_stats']['length']} bp
GC Content: {result['basic_stats']['gc_content']:.2f}%
Molecular Weight: {result['basic_stats']['molecular_weight']:.2f} Da
"""
        text_widget.insert(tk.END, content)

    def check_sequence_quality(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        quality_report = "=== SEQUENCE QUALITY REPORT ===\n\n"
        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            quality_report += f"Sequence: {seq_id}\nLength: {stats['length']} bp\nGC Content: {stats['gc_content']:.2f}%\n"

            issues = []
            if stats['gc_content'] < 20 or stats['gc_content'] > 80:
                issues.append("Extreme GC content")
            if stats['length'] < 100:
                issues.append("Very short sequence")

            quality_report += f"Quality: {'GOOD' if not issues else ', '.join(issues)}\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, quality_report)

    def validate_sequence_format(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        validation_report = "=== SEQUENCE FORMAT VALIDATION ===\n\n"
        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            validation_report += f"Sequence: {seq_id}\n"

            if stats['length'] < 50:
                validation_report += "WARNING: Very short sequence (< 50 bp)\n"
            elif stats['length'] > 100000:
                validation_report += "WARNING: Very long sequence (> 100kb)\n"
            else:
                validation_report += "Length: Valid\n"

            total_bases = sum(stats[base] for base in 'ATGC')
            if total_bases == stats['length']:
                validation_report += "Composition: Pure DNA sequence\n"
            else:
                validation_report += f"WARNING: Contains {stats['length'] - total_bases} non-standard bases\n"

            validation_report += "\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, validation_report)

    def analyze_complexity(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        complexity_report = "=== SEQUENCE COMPLEXITY ANALYSIS ===\n\n"
        for seq_id, result in self.results.items():
            complexity = result.get('complexity', {})
            if complexity:
                complexity_report += f"Sequence: {seq_id}\n"
                complexity_report += f"Shannon Entropy: {complexity['shannon_entropy']:.3f}\n"
                complexity_report += f"Complexity Score: {complexity['complexity_score']:.3f}\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, complexity_report)

    def find_low_complexity(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        low_complexity_report = "=== LOW COMPLEXITY REGIONS ===\n\n"
        for seq_id, result in self.results.items():
            complexity = result.get('complexity', {})
            if complexity and complexity['low_complexity_regions']:
                low_complexity_report += f"Sequence: {seq_id}\n"
                for i, region in enumerate(complexity['low_complexity_regions'][:5], 1):

                    low_complexity_report += f"  Region {i}: {region['start']}-{region['end']} (entropy: {region['entropy']:.3f})\n"
                low_complexity_report += "\n"
            else:
                low_complexity_report += f"Sequence: {seq_id} - No low complexity regions found\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, low_complexity_report)

    def calculate_at_gc_ratios(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        ratio_report = "=== AT/GC RATIO ANALYSIS ===\n\n"
        ratio_report += f"{'Sequence ID':<20} {'AT/GC Ratio':<12} {'Category':<15}\n" + "-" * 50 + "\n"

        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            at_count = stats['A'] + stats['T']
            gc_count = stats['G'] + stats['C']
            at_gc_ratio = at_count / gc_count if gc_count > 0 else float('inf')

            category = "GC-rich" if at_gc_ratio < 0.5 else "AT-rich" if at_gc_ratio > 2.0 else "Balanced"
            ratio_report += f"{seq_id:<20} {at_gc_ratio:<12.2f} {category:<15}\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, ratio_report)

    def sliding_window_analysis(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        window_size = ask_integer("Window Size", "Enter window size (bp):", initial=100, minval=10, maxval=1000)
        if not window_size:
            return

        sliding_report = f"=== SLIDING WINDOW ANALYSIS (Window: {window_size} bp) ===\n\n"
        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            sliding_report += f"Sequence: {seq_id} (Length: {stats['length']} bp)\n"
            sliding_report += f"Average GC content: {stats['gc_content']:.2f}%\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, sliding_report)

    def compare_all_sequences(self):
        if not self.results or len(self.results) < 2:
            messagebox.showwarning("Warning", "Need at least 2 sequences for comparison")
            return

        comparison_report = "=== SEQUENCE COMPARISON REPORT ===\n\n"
        seq_ids = list(self.results.keys())
        comparison_report += f"Comparing {len(seq_ids)} sequences:\n"
        for seq_id in seq_ids:
            comparison_report += f"  - {seq_id}\n"

        lengths = [self.results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]
        gc_contents = [self.results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]

        comparison_report += f"\nLength range: {min(lengths)} - {max(lengths)} bp\n"
        comparison_report += f"Average length: {np.mean(lengths):.0f} bp\n"
        comparison_report += f"GC content range: {min(gc_contents):.2f}% - {max(gc_contents):.2f}%\n"
        comparison_report += f"Average GC content: {np.mean(gc_contents):.2f}%\n"

        self.comparison_text.delete(1.0, tk.END)
        self.comparison_text.insert(tk.END, comparison_report)

    def export_summary_csv(self):
        if not self.results:
            messagebox.showwarning("Warning", "No results to export")
            return

        filename = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if filename:
            try:
                summary_data = []
                for seq_id, result in self.results.items():
                    stats = result['basic_stats']
                    summary_data.append({
                        'Sequence_ID': seq_id,
                        'Length_bp': stats['length'],
                        **{f'{base}_count': stats[base] for base in 'ATGC'},
                        'GC_Content_%': round(stats['gc_content'], 2),
                        'Molecular_Weight_Da': round(stats['molecular_weight'], 2),
                    })

                pd.DataFrame(summary_data).to_csv(filename, index=False)
                messagebox.showinfo("Success", f"Summary exported to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"CSV export failed: {str(e)}")


    def export_detailed_csv(self):
        if not self.results:
            messagebox.showwarning("Warning", "No results to export")
            return

        filename = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if filename:
            try:
                detailed_data = []
                for seq_id, result in self.results.items():
                    stats = result['basic_stats']
                    row = {
                        'Sequence_ID': seq_id,
                        'Length_bp': stats['length'],
                        **{f'{base}_count': stats[base] for base in 'ATGC'},
                        'GC_Content_%': round(stats['gc_content'], 2),
                        'Molecular_Weight_Da': round(stats['molecular_weight'], 2),
                        'Tandem_Repeats': len(result['tandem_repeats']),
                        'Total_Motifs': len(result['motifs'])
                    }

                    if result.get('complexity'):
                        complexity = result['complexity']
                        row.update({
                            'Shannon_Entropy': round(complexity.get('shannon_entropy', 0), 4),
                            'Complexity_Score': round(complexity.get('complexity_score', 0), 4),
                            'Low_Complexity_Regions': len(complexity.get('low_complexity_regions', []))
                        })

                    detailed_data.append(row)

                pd.DataFrame(detailed_data).to_csv(filename, index=False)
                messagebox.showinfo("Success", f"Detailed results exported to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Detailed CSV export failed: {str(e)}")

    def sequence_converter(self):
        sequence_converter_window(self.root)

    def reverse_complement_tool(self):
        sequence = ask_string("Reverse Complement", "Enter DNA sequence:")
        if sequence:
            try:
                result = str(Seq(sequence.upper()).reverse_complement())
                messagebox.showinfo("Result", f"Reverse complement: {result}")
            except Exception as e:
                messagebox.showerror("Error", f"Invalid sequence: {str(e)}")


    def translation_tool(self):
        sequence = ask_string("Translation", "Enter DNA/RNA sequence:")
        if sequence:
            try:
                protein = str(Seq(sequence.upper()).translate())
                messagebox.showinfo("Result", f"Protein: {protein}")
            except Exception as e:
                messagebox.showerror("Error", f"Translation failed: {str(e)}")


    def batch_file_processor(self):
        files = filedialog.askopenfilenames(title='Select multiple FASTA files for batch processing', filetypes=[('FASTA files', '*.fasta *.fa *.fas'), ('All files', '*.*')])
        if not files:
            return

        window = tk.Toplevel(self.root)
        window.title("Batch Processing Results")
        window.geometry("800x600")

        text_widget = tk.Text(window, wrap=tk.WORD)
        text_widget.pack(fill="both", expand=True, padx=10, pady=10)

        batch_report = "=== BATCH PROCESSING SUMMARY ===\n\n"
        total_sequences = 0

        for file_path in files:
            try:
                file_sequences = sum(1 for _ in SeqIO.parse(file_path, 'fasta'))
                total_sequences += file_sequences
                batch_report += f"{os.path.basename(file_path)}: {file_sequences} sequences\n"
            except Exception as e:
                batch_report += f"Error processing {os.path.basename(file_path)}: {str(e)}\n"

        batch_report += f"\nTotal sequences: {total_sequences}\n"
        text_widget.insert(tk.END, batch_report)


    def sequence_length_filter(self):
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        window = tk.Toplevel(self.root)
        window.title("Sequence Length Filter")
        window.geometry("500x400")


        input_frame = tk.Frame(window)
        input_frame.pack(pady=10, padx=20, fill='x')

        tk.Label(input_frame, text="Minimum length (bp):").grid(row=0, column=0, sticky='w', pady=5)
        min_entry = tk.Entry(input_frame, width=10)
        min_entry.insert(0, "0")
        min_entry.grid(row=0, column=1, padx=10)

        tk.Label(input_frame, text="Maximum length (bp):").grid(row=1, column=0, sticky='w', pady=5)
        max_entry = tk.Entry(input_frame, width=10)
        max_entry.insert(0, "999999")
        max_entry.grid(row=1, column=1, padx=10)

        result_text = tk.Text(window, height=15)
        result_text.pack(pady=10, padx=20, fill='both', expand=True)

        def apply_filter():
            try:
                min_length = int(min_entry.get())
                max_length = int(max_entry.get())

                passed = [(seq_id, result['basic_stats']['length']) for seq_id, result in self.results.items()
                          if min_length <= result['basic_stats']['length'] <= max_length]
                failed = [(seq_id, result['basic_stats']['length']) for seq_id, result in self.results.items()
                          if not (min_length <= result['basic_stats']['length'] <= max_length)]

                filter_report = f"=== LENGTH FILTER RESULTS ===\nFilter: {min_length} - {max_length} bp\n\n"
                filter_report += f"PASSED ({len(passed)}):\n" + "\n".join(f"{seq_id}: {length:,} bp" for seq_id, length in passed)
                filter_report += f"\n\nFAILED ({len(failed)}):\n" + "\n".join(f"{seq_id}: {length:,} bp" for seq_id, length in failed)
                filter_report += f"\n\nPass rate: {len(passed)/len(self.results)*100:.1f}%"

                result_text.delete(1.0, tk.END)
                result_text.insert(tk.END, filter_report)
            except ValueError:
                messagebox.showerror("Error", "Please enter valid numbers")


        tk.Button(input_frame, text="Apply Filter", command=apply_filter).grid(row=2, column=0, columnspan=2, pady=10)

    def new_analysis(self):
        self.file_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.results = {}
        self.summary_text.delete(1.0, tk.END)
        self.detailed_text.delete(1.0, tk.END)
        self.progress['value'] = 0
        self.status_label.config(text="Ready")

        for label in self.plot_labels:
            label.destroy()
        self.plot_images = []
        self.plot_labels = []


    def save_session(self):
        filename = filedialog.asksaveasfilename(defaultextension=".session", filetypes=[("Session files", "*.session")])
        if filename:
            session_data = {'files': self.file_entry.get(), 'output': self.output_entry.get(), 'results': self.results}
            save_session(filename, session_data)
            messagebox.showinfo("Success", "Session saved")


    def load_session(self):
        filename = filedialog.askopenfilename(filetypes=[("Session files", "*.session")])
        if filename:
            session_data = load_session(filename)
            if session_data:
                self.file_entry.delete(0, tk.END)
                self.file_entry.insert(0, session_data.get('files', ''))
                self.output_entry.delete(0, tk.END)
                self.output_entry.insert(0, session_data.get('output', ''))
                self.results = session_data.get('results', {})
                if self.results:
                    self.display_results(self.results)
                messagebox.showinfo("Success", "Session loaded")


    def generate_report(self):
        if not self.results:
            messagebox.showwarning("Warning", "No results to generate report")
            return

        output_dir = self.output_entry.get() or filedialog.askdirectory(title="Select directory for report")
        if not output_dir:
            return

        report_type = ask_string("Report Type", "Enter report type (PDF or Excel):", initial="PDF")
        if not report_type:
            return
        report_type = report_type.lower()

        try:
            if report_type == "pdf":
                filename = generate_pdf_report(self.results, output_dir)
            elif report_type == "excel":
                filename = generate_excel_report(self.results, output_dir)
            else:
                messagebox.showwarning("Warning", "Invalid report type. Please choose PDF or Excel.")
                return

            messagebox.showinfo("Success", f"Comprehensive {report_type.upper()} report generated successfully!\nSaved to: {filename}")

        except Exception as e:
            logging.exception("Report generation failed")
            messagebox.showerror("Error", f"Report generation failed: {str(e)}")


    def show_help(self):
        help_text = """Professional RNA Analysis Bioinformatics - User Guide

MAIN FEATURES:
Comprehensive sequence analysis
Quality Control and validation
Sequence comparison tools
Visualizations and reports

USAGE:
1. Select FASTA/GenBank files
2. Choose output directory
3. Click 'Start Analysis'
4. View results in Results tab"""
        messagebox.showinfo("User guide", help_text)


    def show_about(self):
        about_text = """Professional RNA Analysis Bioinformatics

A comprehensive toolkit for DNA/RNA sequence analysis
designed for research and educational purposes.

Built with Python, BioPython, and tkinter"""
        messagebox.showinfo("About", about_text)

if __name__ == "__main__":
    import sys, logging, os
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        # Prefer ttkbootstrap Window so bootstyle is recognized
        root = tb.Window(themename="flatly")
        app = EnhancedDNAToolGUI(root)
        logging.info("GUI initialized; entering mainloop")
        root.mainloop()
    except Exception as e:
        logging.exception("Failed to start GUI")
        # Explicit headless hint for CI/workflow logs
        if os.name != "nt" and not os.environ.get("DISPLAY"):
            print("Error: No display detected ($DISPLAY unset). GUI cannot launch in this environment.")
        sys.exit(1)
