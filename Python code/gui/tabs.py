
"""
Μονάδα: EnhancedDNAToolGUI
Κύρια εφαρμογή για ανάλυση RNA βιοπληροφορικής.
Προσφέρει modular GUI με Tkinter για φόρτωση, ανάλυση, οπτικοποίηση και αναφορές αλληλουχιών DNA/RNA.
"""
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
    """
    Κύρια κλάση που υλοποιεί το γραφικό περιβάλλον χρήστη (GUI) για την ανάλυση βιολογικών αλληλουχιών.
    Διαχειρίζεται την δημιουργία καρτελών (tabs), τη φόρτωση αρχείων, την εκκίνηση αναλύσεων και την εμφάνιση αποτελεσμάτων.
    """

    def __init__(self, root):
        """
        Αρχικοποιεί το GUI, ορίζει τον βασικό παράθυρο και την αρχική κατάσταση.
        
        Παράμετροι:
            root (Tk): Η ρίζα του Tkinter παραθύρου εφαρμογής.
        """
        self.root = root
        root.title('Professional RNA Analysis Bioinformatics')
        root.geometry('1400x900')
        self.session = {}  # Αποθήκευση δεδομένων συνεδρίας
        self.results = {}  # Αποθήκευση αποτελεσμάτων ανάλυσης
        self.setup_widgets()  # Δημιουργία γραφικών στοιχείων

    def setup_widgets(self):
        """
        Δημιουργεί το widget Notebook με τις βασικές καρτέλες της εφαρμογής και τις προσθέτει.
        Στη συνέχεια καλεί συναρτήσεις ρύθμισης για κάθε καρτέλα.
        """
        # Δημιουργία Notebook για καρτέλες
        nb = tb.Notebook(self.root)
        nb.pack(fill='both', expand=True)

        self.frames = {}  # Λεξικό για αποθήκευση των frames κάθε καρτέλας
        
        # Ορισμός καρτελών: (κλειδί, τίτλος)
        tabs = [('main', 'Main Analysis'), ('advanced', 'Advanced Analysis'),
                ('qc', 'Quality Control'), ('compare', 'Sequence Comparison'), ('results', 'Results')]

        # Δημιουργία frame για κάθε καρτέλα
        for key, text in tabs:
            self.frames[key] = tb.Frame(nb)
            nb.add(self.frames[key], text=text)

        # Ρύθμιση περιεχομένου κάθε καρτέλας
        self.setup_main_tab()
        self.setup_advanced_tab()
        self.setup_qc_tab()
        self.setup_comparison_tab()
        self.setup_results_tab()
        self.setup_menu()

    def setup_main_tab(self):
        """
        Δημιουργεί και διατάσσει τα γραφικά στοιχεία (widgets) στην καρτέλα 'Main Analysis'.
        Περιλαμβάνει επιλογή αρχείων FASTA/GenBank, φάκελο εξόδου, progress bar,
        επιλογές δημιουργίας γραφημάτων, και κουμπί έναρξης ολοκληρωμένης ανάλυσης.
        """
        # Ετικέτα για επιλογή αρχείων
        tb.Label(self.frames['main'], text='Select FASTA/GenBank files:', font=('Arial', 12, 'bold')).pack(pady=6)
        
        # Πεδίο εισαγωγής για διαδρομές αρχείων
        self.file_entry = tb.Entry(self.frames['main'], width=120)
        self.file_entry.pack(pady=3)

        # Frame για κουμπιά επιλογής αρχείων
        btn_frame = tb.Frame(self.frames['main'])
        btn_frame.pack(pady=5)
        tb.Button(btn_frame, text='Browse Files', bootstyle='primary', command=self.browse_files).pack(side='left', padx=5)

        # Ετικέτα και πεδίο για φάκελο εξόδου
        tb.Label(self.frames['main'], text='Output folder:', font=('Arial', 12, 'bold')).pack(pady=(10, 6))
        self.output_entry = tb.Entry(self.frames['main'], width=120)
        self.output_entry.pack(pady=3)
        tb.Button(self.frames['main'], text='Select Output Folder', bootstyle='primary', command=self.select_output_folder).pack(pady=3)

        # Progress bar για την πρόοδο ανάλυσης
        self.progress = tb.Progressbar(self.frames['main'], length=600, mode='determinate')
        self.progress.pack(pady=10)
        
        # Ετικέτα κατάστασης
        self.status_label = tb.Label(self.frames['main'], text='Ready', font=('Arial', 10))
        self.status_label.pack(pady=5)

        # Frame για επιλογές γραφημάτων
        toggles_frame = tb.Frame(self.frames['main'])
        toggles_frame.pack(pady=10)

        tb.Label(toggles_frame, text='Plot Generation Options:', font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)

        # Checkbox για ατομικά γραφήματα
        self.individual_plots_var = tk.BooleanVar(value=True)
        individual_check = tb.Checkbutton(toggles_frame, text='Generate Individual Plots', variable=self.individual_plots_var, bootstyle='success')
        individual_check.grid(row=1, column=0, padx=20, pady=5, sticky='w')

        # Checkbox για συνδυαστικά γραφήματα
        self.joint_plots_var = tk.BooleanVar(value=True)
        joint_check = tb.Checkbutton(toggles_frame, text='Generate Joint Plots', variable=self.joint_plots_var, bootstyle='success')
        joint_check.grid(row=1, column=1, padx=20, pady=5, sticky='w')

        # Κουμπί έναρξης ανάλυσης
        tb.Button(self.frames['main'], text='Start Comprehensive Analysis', bootstyle='success-outline', command=self.start_analysis).pack(pady=15)

    def setup_advanced_tab(self):
        """
        Δημιουργεί και διατάσσει τα γραφικά στοιχεία στην καρτέλα 'Advanced Analysis'.
        Περιλαμβάνει επιλογή αλληλουχιών για ευθυγράμμιση (alignment) και πεδίο για
        χειροκίνητη εισαγωγή αλληλουχίας προς ανάλυση.
        """
        # Frame για ευθυγράμμιση αλληλουχιών
        align_frame = tb.LabelFrame(self.frames['advanced'], text='Sequence Alignment', padding=10)
        align_frame.pack(pady=10, padx=20, fill='x')

        tb.Label(align_frame, text='Select two sequences to align:').pack(anchor='w')
        
        # Combobox για επιλογή πρώτης αλληλουχίας
        self.seq1_combo = ttk.Combobox(align_frame, width=50, state='readonly')
        self.seq1_combo.pack(pady=5, fill='x')
        
        # Combobox για επιλογή δεύτερης αλληλουχίας
        self.seq2_combo = ttk.Combobox(align_frame, width=50, state='readonly')
        self.seq2_combo.pack(pady=5, fill='x')
        
        # Κουμπί εκτέλεσης ευθυγράμμισης
        tb.Button(align_frame, text='Perform Alignment', bootstyle='info', command=self.perform_alignment).pack(pady=5)

        # Frame για χειροκίνητη εισαγωγή αλληλουχίας
        manual_frame = tb.LabelFrame(self.frames['advanced'], text='Manual Sequence Input', padding=10)
        manual_frame.pack(pady=10, padx=20, fill='both', expand=True)

        tb.Label(manual_frame, text='Paste sequence:').pack(anchor='w')
        
        # Πεδίο κειμένου για εισαγωγή αλληλουχίας
        self.manual_seq_text = tk.Text(manual_frame, height=8, width=80)
        self.manual_seq_text.pack(pady=5, fill='both', expand=True)

        # Frame για κουμπιά χειροκίνητης ανάλυσης
        manual_btn_frame = tb.Frame(manual_frame)
        manual_btn_frame.pack(pady=5)

        # Κουμπιά για ανάλυση και καθαρισμό
        tb.Button(manual_btn_frame, text='Analyze Single Sequence', bootstyle='warning', command=self.analyze_manual_sequence).pack(side='left', padx=5)
        tb.Button(manual_btn_frame, text='Clear', bootstyle='secondary', command=lambda: self.manual_seq_text.delete(1.0, tk.END)).pack(side='left', padx=5)

    def setup_qc_tab(self):
        """
        Ρυθμίζει την καρτέλα 'Quality Control' με τα κουμπιά για λειτουργίες
        ελέγχου ποιότητας αλληλουχιών και παράθεση των αποθηκευμένων πληροφοριών
        σε πεδίο κειμένου.
        """
        # Τίτλος καρτέλας
        tb.Label(self.frames['qc'], text='Sequence Quality Control', font=('Arial', 14, 'bold')).pack(pady=10)

        # Λίστα κουμπιών ελέγχου ποιότητας: (κείμενο, στυλ, εντολή)
        qc_buttons = [
            ('Check Sequence Quality', 'primary', self.check_sequence_quality),
            ('Validate Sequence Format', 'warning', self.validate_sequence_format),
            ('Calculate AT/GC Ratio', 'info', self.calculate_at_gc_ratios),
            ('Sliding Window Analysis', 'success', self.sliding_window_analysis),
            ('Analyze Sequence Complexity', 'info', self.analyze_complexity),
            ('Find Low Complexity Regions', 'secondary', self.find_low_complexity)
        ]

        # Δημιουργία κουμπιών
        for text, style, command in qc_buttons:
            tb.Button(self.frames['qc'], text=text, bootstyle=style, command=command).pack(pady=2)

        # Πεδίο κειμένου για αποτελέσματα ελέγχου ποιότητας
        self.qc_text = tk.Text(self.frames['qc'], height=15, width=100)
        self.qc_text.pack(pady=10, padx=20, fill='both', expand=True)

    def setup_comparison_tab(self):
        """
        Ρυθμίζει την καρτέλα 'Sequence Comparison' με εργαλεία για σύγκριση
        ακολουθιών και εμφάνιση αποτελεσμάτων.
        """
        # Τίτλος καρτέλας
        tb.Label(self.frames['compare'], text='Sequence Comparison Tools', font=('Arial', 14, 'bold')).pack(pady=10)

        # Frame για πολλαπλή ευθυγράμμιση αλληλουχιών
        msa_frame = tb.LabelFrame(self.frames['compare'], text='Multiple Sequence Alignment', padding=10)
        msa_frame.pack(pady=10, padx=20, fill='x')
        tb.Button(msa_frame, text='Compare All Sequences', bootstyle='info', command=self.compare_all_sequences).pack(pady=5)

        # Πεδίο κειμένου για αποτελέσματα σύγκρισης
        self.comparison_text = tk.Text(self.frames['compare'], height=20, width=100)
        self.comparison_text.pack(pady=10, padx=20, fill='both', expand=True)

    def setup_results_tab(self):
        """
        Δημιουργεί την καρτέλα 'Results' με υπο-καρτέλες για συνοπτική παρουσίαση,
        λεπτομερή αποτελέσματα και οπτικοποιήσεις.
        """
        # Notebook για υπο-καρτέλες αποτελεσμάτων
        results_nb = tb.Notebook(self.frames['results'])
        results_nb.pack(fill='both', expand=True, padx=10, pady=10)

        # Υπο-καρτέλα συνοπτικών αποτελεσμάτων
        self.summary_frame = tb.Frame(results_nb)
        results_nb.add(self.summary_frame, text='Summary')

        # Κουμπί εξαγωγής συνοπτικών σε CSV
        summary_button_frame = tb.Frame(self.summary_frame)
        summary_button_frame.pack(pady=5)
        tb.Button(summary_button_frame, text='Export Summary to CSV', bootstyle='info', command=self.export_summary_csv).pack(side='left', padx=5)

        # Πεδίο κειμένου για συνοπτικά αποτελέσματα
        self.summary_text = tk.Text(self.summary_frame, height=25, width=100)
        self.summary_text.pack(fill='both', expand=True, padx=5, pady=5)

        # Υπο-καρτέλα λεπτομερών αποτελεσμάτων
        self.detailed_frame = tb.Frame(results_nb)
        results_nb.add(self.detailed_frame, text='Detailed Results')

        # Κουμπί εξαγωγής λεπτομερών σε CSV
        detailed_button_frame = tb.Frame(self.detailed_frame)
        detailed_button_frame.pack(pady=5)
        tb.Button(detailed_button_frame, text='Export Detailed Results to CSV', bootstyle='info', command=self.export_detailed_csv).pack(side='left', padx=5)

        # Πεδίο κειμένου για λεπτομερή αποτελέσματα
        self.detailed_text = tk.Text(self.detailed_frame, height=25, width=100)
        self.detailed_text.pack(fill='both', expand=True, padx=5, pady=5)

        # Υπο-καρτέλα οπτικοποιήσεων
        self.plots_frame = tb.Frame(results_nb)
        results_nb.add(self.plots_frame, text='Visualizations')

        # Ρύθμιση στοιχείων ελέγχου οπτικοποίησης
        self.setup_visualization_controls()

        # Λίστες για αποθήκευση εικόνων και ετικετών γραφημάτων
        self.plot_images = []
        self.plot_labels = []

    def setup_menu(self):
        """
        Δημιουργεί το κύριο μενού της εφαρμογής με κατηγορίες 'File', 'Tools' και 'Help'
        που περιλαμβάνουν διάφορες ενέργειες.
        """
        menubar = tk.Menu(self.root)

        # Μενού File
        file_menu = tk.Menu(menubar, tearoff=0)
        for item in [('New Analysis', self.new_analysis), ('Save Session', self.save_session),
                     ('Load Session', self.load_session), ('---', None), ('Exit', self.root.quit)]:
            if item[1]:
                file_menu.add_command(label=item[0], command=item[1])
            else:
                file_menu.add_separator()
        menubar.add_cascade(label='File', menu=file_menu)

        # Μενού Tools
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

        # Μενού Help
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label='User Guide', command=self.show_help)
        help_menu.add_command(label='About', command=self.show_about)
        menubar.add_cascade(label='Help', menu=help_menu)

        self.root.config(menu=menubar)

    def browse_files(self):
        """
        Άνοιγμα διαλόγου για επιλογή αρχείων FASTA ή GenBank.
        Τα επιλεγμένα αρχεία προστίθενται στην είσοδο αρχείων.
        """
        files = filedialog.askopenfilenames(
            title='Select FASTA/GenBank files',
            filetypes=[('FASTA files', '*.fasta *.fa *.fas'), ('GenBank files', '*.gb *.gbk'), ('All supported', '*.fasta *.fa *.fas *.gb *.gbk')]
        )
        if files:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, ';'.join(files))

    def select_output_folder(self):
        """
        Άνοιγμα διαλόγου για επιλογή φακέλου εξόδου.
        Ο επιλεγμένος φάκελος προστίθεται στο πεδίο φακέλου εξόδου.
        """
        folder = filedialog.askdirectory(title='Select output folder')
        if folder:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, folder)

    def start_analysis(self):
        """
        Εκκίνηση της διαδικασίας ανάλυσης σε ξεχωριστό νήμα.
        Ελέγχει αν έχουν επιλεγεί αρχεία και φάκελος εξόδου πριν ξεκινήσει.
        """
        # Λήψη επιλεγμένων αρχείων και φακέλου εξόδου
        files = self.file_entry.get().split(';') if self.file_entry.get() else []
        output_dir = self.output_entry.get()

        # Έλεγχος αν υπάρχουν αρχεία και φάκελος εξόδου
        if not files or not output_dir:
            messagebox.showwarning('Warning', 'Please select files and output directory')
            return

        # Δημιουργία και εκκίνηση νήματος ανάλυσης
        thread = threading.Thread(target=self.run_analysis_thread, args=(files, output_dir))
        thread.daemon = True
        thread.start()

    def run_analysis_thread(self, files, output_dir):
        """
        Εκτελεί την ανάλυση σε ξεχωριστό νήμα για να μην παγώνει το GUI.
        Επεξεργάζεται όλα τα αρχεία, εκτελεί ανάλυση και δημιουργεί αναφορές.
        
        Παράμετροι:
            files (list): Λίστα με διαδρομές αρχείων προς ανάλυση
            output_dir (str): Φάκελος εξόδου για αποτελέσματα
        """
        try:
            import gc
            start_time = time.time()  # Έναρξη χρονομέτρησης
            self.update_status("Starting analysis...")

            # Έλεγχος λειτουργίας χαμηλής μνήμης
            is_low_memory = low_memory_mode()
            if is_low_memory:
                self.update_status("Low memory mode enabled - some features disabled")

            # Μέτρηση συνολικών αλληλουχιών
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

            # Αρχικοποίηση μεταβλητών για αποτελέσματα
            all_results = {}
            seq_names = []
            total_processed = 0

            # Επεξεργασία κάθε αλληλουχίας
            for seq_id, seq in self.parse_files(files):
                total_processed += 1
                self.update_status(f"Analyzing {seq_id} ({total_processed}/{total_sequences})...")

                # Παράλειψη μεγάλων αλληλουχιών σε λειτουργία χαμηλής μνήμης
                if is_low_memory and len(seq) > 50000:
                    logging.warning(f"Skipping large sequence {seq_id} in low memory mode")
                    self.update_progress((total_processed / total_sequences) * 80)
                    continue

                try:
                    # Εκτέλεση ολοκληρωμένης ανάλυσης
                    result = comprehensive_analysis(seq, seq_id)
                    all_results[seq_id] = result
                    seq_names.append(seq_id)

                    # Απελευθέρωση μνήμης
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

            # Ενημέρωση combobox με ονόματα αλληλουχιών
            self.root.after(0, lambda: self.update_sequence_combos(seq_names))
            self.results = all_results

            # Δημιουργία αναφορών και οπτικοποιήσεων
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

            # Οριστικοποίηση αποτελεσμάτων
            self.update_status("Finalizing results...")
            self.update_progress(95)
            self.root.after(0, lambda: self.display_results(all_results))
            self.root.after(0, lambda: self.load_plots())
            self.update_progress(100)

            # Τέλος χρονομέτρησης
            end_time = time.time()
            elapsed = end_time - start_time
            elapsed_str = f"{elapsed:.2f} seconds ({elapsed/60:.2f} minutes)"

            # Μόνιμη καταγραφή
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

            # Ενημέρωση GUI
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
        """
        Αναλύει αρχεία FASTA/GenBank και επιστρέφει generator με (id, sequence).
        
        Παράμετροι:
            files (list): Λίστα διαδρομών αρχείων
            
        Yields:
            tuple: (seq_id, seq) για κάθε αλληλουχία
        """
        for file_path in files:
            try:
                # Προσδιορισμός τύπου αρχείου
                format_type = 'genbank' if file_path.lower().endswith(('.gb', '.gbk')) else 'fasta'
                for record in SeqIO.parse(file_path, format_type):
                    yield (record.id, record.seq)
            except Exception as e:
                logging.error(f"Failed to parse {file_path}: {e}")

    def update_status(self, message):
        """
        Ενημερώνει την ετικέτα κατάστασης με το δεδομένο μήνυμα.
        
        Παράμετροι:
            message (str): Μήνυμα κατάστασης
        """
        def update():
            self.status_label.config(text=message)
            logging.info(message)
        self.root.after(0, update)

    def update_progress(self, value):
        """
        Ενημερώνει την μπάρα προόδου με τη δεδομένη τιμή.
        
        Παράμετροι:
            value (float): Τιμή προόδου (0-100)
        """
        def update():
            self.progress['value'] = value
        self.root.after(0, update)

    def update_sequence_combos(self, seq_names):
        """
        Ενημερώνει τα combobox επιλογής αλληλουχιών με τα ονόματα αλληλουχιών.
        
        Παράμετροι:
            seq_names (list): Λίστα με ονόματα αλληλουχιών
        """
        self.seq1_combo['values'] = seq_names
        self.seq2_combo['values'] = seq_names

    def display_results(self, results):
        """
        Εμφανίζει τα αποτελέσματα ανάλυσης στα πεδία κειμένου.
        
        Παράμετροι:
            results (dict): Λεξικό με αποτελέσματα ανάλυσης
        """
        # Καθαρισμός πεδίων κειμένου
        self.summary_text.delete(1.0, tk.END)
        self.detailed_text.delete(1.0, tk.END)

        # Δημιουργία συνοπτικής αναφοράς
        summary = "=== ANALYSIS SUMMARY ===\n\n"
        for seq_id, result in results.items():
            stats = result['basic_stats']
            summary += f"Sequence: {seq_id}\n  Length: {stats['length']} bp\n  GC Content: {stats['gc_content']:.2f}%\n\n"
        self.summary_text.insert(tk.END, summary)

        # Δημιουργία λεπτομερούς αναφοράς
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
        """
        Ρυθμίζει τα στοιχεία ελέγχου για την επιλογή τρόπου οπτικοποίησης
        (συνδυαστικά ή ατομικά γραφήματα).
        """
        # Frame για στοιχεία ελέγχου
        controls_frame = tb.Frame(self.plots_frame)
        controls_frame.pack(fill='x', padx=10, pady=5)

        # Επιλογή τρόπου οπτικοποίησης
        tb.Label(controls_frame, text='Visualization Mode:', font=('Arial', 12, 'bold')).pack(side='left', padx=(0, 10))

        self.viz_mode_var = tk.StringVar(value='Joint Plots')
        self.viz_mode_combo = tb.Combobox(controls_frame, textvariable=self.viz_mode_var,
                                          values=['Joint Plots', 'Individual Plots'],
                                          state='readonly', width=20)
        self.viz_mode_combo.pack(side='left', padx=(0, 10))
        self.viz_mode_combo.bind('<<ComboboxSelected>>', self.on_viz_mode_change)

        # Frame για επιλογή φακέλου ατομικών γραφημάτων (αρχικά κρυφό)
        self.folder_selection_frame = tb.Frame(controls_frame)
        self.folder_selection_frame.pack(side='left', padx=(20, 0))

        self.folder_label = tb.Label(self.folder_selection_frame, text='Select Sequence:', font=('Arial', 10, 'bold'))
        self.folder_label.pack(side='left', padx=(0, 5))

        self.folder_var = tk.StringVar()
        self.folder_combo = tb.Combobox(self.folder_selection_frame, textvariable=self.folder_var,
                                        state='readonly', width=25)
        self.folder_combo.pack(side='left', padx=(0, 10))
        self.folder_combo.bind('<<ComboboxSelected>>', self.on_folder_change)

        # Κουμπιά πλοήγησης για ατομικά γραφήματα
        self.nav_frame = tb.Frame(self.folder_selection_frame)
        self.nav_frame.pack(side='left', padx=(10, 0))

        self.prev_button = tb.Button(self.nav_frame, text='← Previous', bootstyle='secondary',
                                     command=self.previous_folder, width=10)
        self.prev_button.pack(side='left', padx=(0, 5))

        self.next_button = tb.Button(self.nav_frame, text='Next →', bootstyle='secondary',
                                     command=self.next_folder, width=10)
        self.next_button.pack(side='left')

        # Αρχικά κρύβει τα στοιχεία ελέγχου ατομικών γραφημάτων
        self.folder_selection_frame.pack_forget()

        # Scrollable frame για γραφήματα
        self.plots_canvas = tb.Canvas(self.plots_frame)
        self.plots_scrollbar = tb.Scrollbar(self.plots_frame, orient='vertical', command=self.plots_canvas.yview)
        self.scrollable_plots_frame = tb.Frame(self.plots_canvas)

        self.scrollable_plots_frame.bind('<Configure>', lambda e: self.plots_canvas.configure(scrollregion=self.plots_canvas.bbox('all')))

        self.plots_canvas.create_window((0, 0), window=self.scrollable_plots_frame, anchor='nw')
        self.plots_canvas.configure(yscrollcommand=self.plots_scrollbar.set)

        self.plots_canvas.pack(side='left', fill='both', expand=True, padx=(10, 0), pady=10)
        self.plots_scrollbar.pack(side='right', fill='y', pady=10)

    def on_viz_mode_change(self, event=None):
        """
        Χειρίζεται την αλλαγή τρόπου οπτικοποίησης.
        Εμφανίζει ή κρύβει τα στοιχεία ελέγχου ατομικών γραφημάτων.
        """
        mode = self.viz_mode_var.get()
        if mode == 'Individual Plots':
            self.folder_selection_frame.pack(side='left', padx=(20, 0))
            self.update_folder_list()
        else:
            self.folder_selection_frame.pack_forget()
        self.load_plots()

    def update_folder_list(self):
        """
        Ενημερώνει τη λίστα διαθέσιμων φακέλων ατομικών γραφημάτων.
        """
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
        """
        Χειρίζεται την αλλαγή επιλογής φακέλου ατομικών γραφημάτων.
        """
        self.load_plots()

    def previous_folder(self):
        """
        Πλοηγείται στον προηγούμενο φάκελο ατομικών γραφημάτων.
        """
        current_folders = list(self.folder_combo['values'])
        current_folder = self.folder_var.get()
        if current_folder in current_folders:
            current_index = current_folders.index(current_folder)
            if current_index > 0:
                self.folder_var.set(current_folders[current_index - 1])
                self.load_plots()

    def next_folder(self):
        """
        Πλοηγείται στον επόμενο φάκελο ατομικών γραφημάτων.
        """
        current_folders = list(self.folder_combo['values'])
        current_folder = self.folder_var.get()
        if current_folder in current_folders:
            current_index = current_folders.index(current_folder)
            if current_index < len(current_folders) - 1:
                self.folder_var.set(current_folders[current_index + 1])
                self.load_plots()

    def load_plots(self):
        """
        Φορτώνει γραφήματα με βάση τον επιλεγμένο τρόπο οπτικοποίησης
        (συνδυαστικά ή ατομικά).
        """
        # Καθαρισμός υπαρχόντων γραφημάτων
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
        """
        Φορτώνει και εμφανίζει συνδυαστικά γραφήματα.
        
        Παράμετροι:
            output_dir (str): Φάκελος εξόδου όπου βρίσκονται τα γραφήματα
        """
        joint_plots_dir = os.path.join(output_dir, 'plots', 'joint')
        if not os.path.exists(joint_plots_dir):
            # Εμφάνιση μηνύματος αν δεν υπάρχουν συνδυαστικά γραφήματα
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text='No joint plots available. Make sure "Generate Joint Plots" is enabled and run analysis.',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        row, col = 0, 0
        max_cols = 2  # Εμφάνιση 2 γραφημάτων ανά σειρά για καλύτερη σαφήνεια

        # Γραφήματα προτεραιότητας για καλύτερη οργάνωση
        priority_plots = ['dashboard.png', 'gc_comparison.png', 'length_distribution.png',
                          'composition_heatmap.png', 'nucleotide_boxplot.png', 'molecular_weight.png',
                          'correlation_matrix.png']

        # Φόρτωση γραφημάτων προτεραιότητας πρώτα
        for plot_file in priority_plots:
            image_path = os.path.join(joint_plots_dir, plot_file)
            if os.path.exists(image_path):
                self.load_and_display_plot(image_path, plot_file, row, col, 600, 450)
                col += 1
                if col >= max_cols:
                    col = 0
                    row += 1

        # Φόρτωση υπολοίπων γραφημάτων
        for plot_file in sorted(os.listdir(joint_plots_dir)):
            if plot_file.endswith('.png') and plot_file not in priority_plots:
                image_path = os.path.join(joint_plots_dir, plot_file)
                self.load_and_display_plot(image_path, plot_file, row, col, 600, 450)
                col += 1
                if col >= max_cols:
                    col = 0
                    row += 1

    def load_individual_plots(self, output_dir):
        """
        Φορτώνει και εμφανίζει ατομικά γραφήματα για επιλεγμένη αλληλουχία.
        
        Παράμετροι:
            output_dir (str): Φάκελος εξόδου όπου βρίσκονται τα γραφήματα
        """
        individual_plots_dir = os.path.join(output_dir, 'plots', 'individual')
        if not os.path.exists(individual_plots_dir):
            # Εμφάνιση μηνύματος αν δεν υπάρχουν ατομικά γραφήματα
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text='No individual plots available. Make sure "Generate Individual Plots" is enabled and run analysis.',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        selected_folder = self.folder_var.get()
        if not selected_folder:
            # Εμφάνιση μηνύματος για επιλογή αλληλουχίας
            select_msg_label = tb.Label(self.scrollable_plots_frame,
                                        text='Please select a sequence from the dropdown above to view its individual plots.',
                                        font=('Arial', 12), foreground='blue')
            select_msg_label.pack(pady=20)
            self.plot_labels.append(select_msg_label)
            return

        seq_folder_path = os.path.join(individual_plots_dir, selected_folder)
        if not os.path.exists(seq_folder_path):
            # Εμφάνιση μηνύματος σφάλματος
            error_label = tb.Label(self.scrollable_plots_frame,
                                   text=f'Plots for sequence "{selected_folder}" not found.',
                                   font=('Arial', 12), foreground='red')
            error_label.pack(pady=20)
            self.plot_labels.append(error_label)
            return

        # Εμφάνιση κεφαλίδας με πληροφορίες αλληλουχίας
        header_label = tb.Label(self.scrollable_plots_frame,
                                text=f'Individual Analysis Plots for: {selected_folder}',
                                font=('Arial', 14, 'bold'), foreground='navy')
        header_label.pack(pady=(10, 20))
        self.plot_labels.append(header_label)

        # Φόρτωση γραφημάτων για την επιλεγμένη αλληλουχία
        plot_files = [f for f in os.listdir(seq_folder_path) if f.endswith('.png')]
        if not plot_files:
            no_plots_label = tb.Label(self.scrollable_plots_frame,
                                      text=f'No plot files found for sequence "{selected_folder}".',
                                      font=('Arial', 12), foreground='gray')
            no_plots_label.pack(pady=20)
            self.plot_labels.append(no_plots_label)
            return

        # Εμφάνιση γραφημάτων σε μία στήλη για ατομική αλληλουχία (ευκολότερη προβολή)
        for i, plot_file in enumerate(sorted(plot_files)):
            image_path = os.path.join(seq_folder_path, plot_file)
            self.load_and_display_plot(image_path, plot_file, i, 0, 800, 600, single_column=True)

    def load_and_display_plot(self, image_path, plot_name, row, col, width=400, height=300, single_column=False):
        """
        Φορτώνει και εμφανίζει μια εικόνα γραφήματος.
        
        Παράμετροι:
            image_path (str): Διαδρομή του αρχείου εικόνας
            plot_name (str): Όνομα γραφήματος
            row (int): Σειρά στο grid
            col (int): Στήλη στο grid
            width (int): Πλάτος εικόνας
            height (int): Ύψος εικόνας
            single_column (bool): Αν True, εμφάνιση σε μία στήλη
        """
        try:
            # Φόρτωση και αλλαγή μεγέθους εικόνας
            img = Image.open(image_path).resize((width, height))
            photo = ImageTk.PhotoImage(img)
            self.plot_images.append(photo)

            # Δημιουργία frame για γράφημα με τίτλο
            plot_frame = tb.Frame(self.scrollable_plots_frame)
            if single_column:
                plot_frame.pack(pady=10, padx=20, fill='x')
            else:
                plot_frame.grid(row=row, column=col, padx=10, pady=10, sticky='nsew')

            # Προσθήκη τίτλου γραφήματος
            title = plot_name.replace('.png', '').replace('_', ' ').title()
            title_label = tb.Label(plot_frame, text=title, font=('Arial', 11, 'bold'))
            title_label.pack(pady=(0, 5))

            # Προσθήκη εικόνας γραφήματος
            img_label = tb.Label(plot_frame, image=photo)
            img_label.pack()

            self.plot_labels.extend([plot_frame, title_label, img_label])

        except Exception as e:
            logging.error(f"Error loading plot {plot_name}: {e}")

    def perform_alignment(self):
        """
        Εκτελεί ευθυγράμμιση μεταξύ δύο επιλεγμένων αλληλουχιών.
        Ζητά από τον χρήστη να επιλέξει μέθοδο (Needleman-Wunsch ή Smith-Waterman).
        """
        seq1_name = self.seq1_combo.get()
        seq2_name = self.seq2_combo.get()

        # Έλεγχος έγκυρης επιλογής
        if not seq1_name or not seq2_name or seq1_name == seq2_name:
            messagebox.showwarning("Warning", "Please select two different sequences")
            return

        from analysis.alignments import sequence_alignment

        # Ερώτηση χρήστη για μέθοδο
        method = ask_string("Alignment Method", "Choose method (NW = Needleman-Wunsch, SW = Smith-Waterman):", initial="NW")
        if not method:
            return
        method = method.strip().lower()
        if method not in ("nw", "sw"):
            messagebox.showerror("Error", "Invalid method. Please enter NW or SW.")
            return

        # Ανάκτηση αλληλουχιών από αρχικά αρχεία
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

        # Εκτέλεση ευθυγράμμισης
        alignment_result = sequence_alignment(seq1, seq2, method=method)

        if not alignment_result:
            messagebox.showerror("Error", "Sequence alignment failed.")
            return

        # Εμφάνιση αποτελεσμάτων
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
        """
        Αναλύει χειροκίνητα εισαγόμενη αλληλουχία DNA/RNA.
        Καθαρίζει την είσοδο και εκτελεί ολοκληρωμένη ανάλυση.
        """
        sequence_text = self.manual_seq_text.get(1.0, tk.END).strip()
        if not sequence_text:
            messagebox.showwarning("Warning", "Please enter a sequence")
            return

        # Καθαρισμός αλληλουχίας - κρατά μόνο έγκυρες βάσεις
        sequence = re.sub(r'[^ATGCN]', '', sequence_text.upper())
        if not sequence:
            messagebox.showwarning("Warning", "No valid DNA/RNA sequence found")
            return

        # Εκτέλεση ανάλυσης
        seq_obj = Seq(sequence)
        result = comprehensive_analysis(seq_obj, "Manual_Input")

        # Εμφάνιση αποτελεσμάτων σε νέο παράθυρο
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
        """
        Ελέγχει την ποιότητα των φορτωμένων αλληλουχιών.
        Αναφέρει προβλήματα όπως ακραίο GC περιεχόμενο ή πολύ μικρό μήκος.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        quality_report = "=== SEQUENCE QUALITY REPORT ===\n\n"
        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            quality_report += f"Sequence: {seq_id}\nLength: {stats['length']} bp\nGC Content: {stats['gc_content']:.2f}%\n"

            # Εντοπισμός προβλημάτων ποιότητας
            issues = []
            if stats['gc_content'] < 20 or stats['gc_content'] > 80:
                issues.append("Extreme GC content")
            if stats['length'] < 100:
                issues.append("Very short sequence")

            quality_report += f"Quality: {'GOOD' if not issues else ', '.join(issues)}\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, quality_report)

    def validate_sequence_format(self):
        """
        Επικυρώνει τη μορφή των φορτωμένων αλληλουχιών.
        Ελέγχει για ακραία μήκη και μη-τυπικές βάσεις.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        validation_report = "=== SEQUENCE FORMAT VALIDATION ===\n\n"
        for seq_id, result in self.results.items():
            stats = result['basic_stats']
            validation_report += f"Sequence: {seq_id}\n"

            # Έλεγχος μήκους
            if stats['length'] < 50:
                validation_report += "WARNING: Very short sequence (< 50 bp)\n"
            elif stats['length'] > 100000:
                validation_report += "WARNING: Very long sequence (> 100kb)\n"
            else:
                validation_report += "Length: Valid\n"

            # Έλεγχος σύνθεσης
            total_bases = sum(stats[base] for base in 'ATGC')
            if total_bases == stats['length']:
                validation_report += "Composition: Pure DNA sequence\n"
            else:
                validation_report += f"WARNING: Contains {stats['length'] - total_bases} non-standard bases\n"

            validation_report += "\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, validation_report)

    def analyze_complexity(self):
        """
        Αναλύει την πολυπλοκότητα των φορτωμένων αλληλουχιών.
        Υπολογίζει Shannon entropy και complexity score.
        """
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
        """
        Εντοπίζει περιοχές χαμηλής πολυπλοκότητας στις φορτωμένες αλληλουχίες.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        low_complexity_report = "=== LOW COMPLEXITY REGIONS ===\n\n"
        for seq_id, result in self.results.items():
            complexity = result.get('complexity', {})
            if complexity and complexity['low_complexity_regions']:
                low_complexity_report += f"Sequence: {seq_id}\n"
                # Εμφάνιση μέχρι 5 περιοχών
                for i, region in enumerate(complexity['low_complexity_regions'][:5], 1):
                    low_complexity_report += f"  Region {i}: {region['start']}-{region['end']} (entropy: {region['entropy']:.3f})\n"
                low_complexity_report += "\n"
            else:
                low_complexity_report += f"Sequence: {seq_id} - No low complexity regions found\n\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, low_complexity_report)

    def calculate_at_gc_ratios(self):
        """
        Υπολογίζει τους λόγους AT/GC για όλες τις φορτωμένες αλληλουχίες.
        Κατηγοριοποιεί τις αλληλουχίες ως GC-rich, AT-rich ή Balanced.
        """
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

            # Κατηγοριοποίηση
            category = "GC-rich" if at_gc_ratio < 0.5 else "AT-rich" if at_gc_ratio > 2.0 else "Balanced"
            ratio_report += f"{seq_id:<20} {at_gc_ratio:<12.2f} {category:<15}\n"

        self.qc_text.delete(1.0, tk.END)
        self.qc_text.insert(tk.END, ratio_report)

    def sliding_window_analysis(self):
        """
        Εκτελεί ανάλυση ολισθαίνοντος παραθύρου στις φορτωμένες αλληλουχίες.
        Ζητά από τον χρήστη το μέγεθος του παραθύρου.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Ερώτηση για μέγεθος παραθύρου
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
        """
        Συγκρίνει όλες τις φορτωμένες αλληλουχίες.
        Παρέχει στατιστικά για μήκος και GC περιεχόμενο.
        """
        if not self.results or len(self.results) < 2:
            messagebox.showwarning("Warning", "Need at least 2 sequences for comparison")
            return

        comparison_report = "=== SEQUENCE COMPARISON REPORT ===\n\n"
        seq_ids = list(self.results.keys())
        comparison_report += f"Comparing {len(seq_ids)} sequences:\n"
        for seq_id in seq_ids:
            comparison_report += f"  - {seq_id}\n"

        # Υπολογισμός στατιστικών
        lengths = [self.results[seq_id]['basic_stats']['length'] for seq_id in seq_ids]
        gc_contents = [self.results[seq_id]['basic_stats']['gc_content'] for seq_id in seq_ids]

        comparison_report += f"\nLength range: {min(lengths)} - {max(lengths)} bp\n"
        comparison_report += f"Average length: {np.mean(lengths):.0f} bp\n"
        comparison_report += f"GC content range: {min(gc_contents):.2f}% - {max(gc_contents):.2f}%\n"
        comparison_report += f"Average GC content: {np.mean(gc_contents):.2f}%\n"

        self.comparison_text.delete(1.0, tk.END)
        self.comparison_text.insert(tk.END, comparison_report)

    def export_summary_csv(self):
        """
        Εξάγει συνοπτικά αποτελέσματα σε αρχείο CSV.
        Περιλαμβάνει βασικά στατιστικά για κάθε αλληλουχία.
        """
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
        """
        Εξάγει λεπτομερή αποτελέσματα σε αρχείο CSV.
        Περιλαμβάνει πλήρη ανάλυση για κάθε αλληλουχία.
        """
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

                    # Προσθήκη δεδομένων πολυπλοκότητας αν υπάρχουν
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
        """
        Ανοίγει το παράθυρο μετατροπής αλληλουχιών.
        """
        sequence_converter_window(self.root)

    def reverse_complement_tool(self):
        """
        Εργαλείο για υπολογισμό αντίστροφης συμπληρωματικής αλληλουχίας.
        Ζητά αλληλουχία DNA από τον χρήστη και εμφανίζει το αποτέλεσμα.
        """
        sequence = ask_string("Reverse Complement", "Enter DNA sequence:")
        if sequence:
            try:
                result = str(Seq(sequence.upper()).reverse_complement())
                messagebox.showinfo("Result", f"Reverse complement: {result}")
            except Exception as e:
                messagebox.showerror("Error", f"Invalid sequence: {str(e)}")

    def translation_tool(self):
        """
        Εργαλείο για μετάφραση αλληλουχίας DNA/RNA σε πρωτεΐνη.
        Ζητά αλληλουχία από τον χρήστη και εμφανίζει την πρωτεϊνική αλληλουχία.
        """
        sequence = ask_string("Translation", "Enter DNA/RNA sequence:")
        if sequence:
            try:
                protein = str(Seq(sequence.upper()).translate())
                messagebox.showinfo("Result", f"Protein: {protein}")
            except Exception as e:
                messagebox.showerror("Error", f"Translation failed: {str(e)}")

    def batch_file_processor(self):
        """
        Επεξεργάζεται πολλαπλά αρχεία FASTA ταυτόχρονα.
        Εμφανίζει συνοπτική αναφορά για τον αριθμό αλληλουχιών σε κάθε αρχείο.
        """
        files = filedialog.askopenfilenames(title='Select multiple FASTA files for batch processing', filetypes=[('FASTA files', '*.fasta *.fa *.fas'), ('All files', '*.*')])
        if not files:
            return

        # Δημιουργία παραθύρου αποτελεσμάτων
        window = tk.Toplevel(self.root)
        window.title("Batch Processing Results")
        window.geometry("800x600")

        text_widget = tk.Text(window, wrap=tk.WORD)
        text_widget.pack(fill="both", expand=True, padx=10, pady=10)

        batch_report = "=== BATCH PROCESSING SUMMARY ===\n\n"
        total_sequences = 0

        # Επεξεργασία κάθε αρχείου
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
        """
        Φιλτράρει αλληλουχίες με βάση το μήκος τους.
        Επιτρέπει στον χρήστη να ορίσει ελάχιστο και μέγιστο μήκος.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No sequences loaded")
            return

        # Δημιουργία παραθύρου φίλτρου
        window = tk.Toplevel(self.root)
        window.title("Sequence Length Filter")
        window.geometry("500x400")

        # Frame για πεδία εισαγωγής
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

        # Πεδίο κειμένου για αποτελέσματα
        result_text = tk.Text(window, height=15)
        result_text.pack(pady=10, padx=20, fill='both', expand=True)

        def apply_filter():
            """
            Εφαρμόζει το φίλτρο μήκους και εμφανίζει τα αποτελέσματα.
            """
            try:
                min_length = int(min_entry.get())
                max_length = int(max_entry.get())

                # Διαχωρισμός αλληλουχιών που περνούν και αποτυγχάνουν το φίλτρο
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

        # Κουμπί εφαρμογής φίλτρου
        tk.Button(input_frame, text="Apply Filter", command=apply_filter).grid(row=2, column=0, columnspan=2, pady=10)

    def new_analysis(self):
        """
        Ξεκινά νέα ανάλυση καθαρίζοντας όλα τα δεδομένα της τρέχουσας συνεδρίας.
        """
        self.file_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.results = {}
        self.summary_text.delete(1.0, tk.END)
        self.detailed_text.delete(1.0, tk.END)
        self.progress['value'] = 0
        self.status_label.config(text="Ready")

        # Καθαρισμός γραφημάτων
        for label in self.plot_labels:
            label.destroy()
        self.plot_images = []
        self.plot_labels = []

    def save_session(self):
        """
        Αποθηκεύει την τρέχουσα συνεδρία σε αρχείο.
        Περιλαμβάνει επιλεγμένα αρχεία, φάκελο εξόδου και αποτελέσματα.
        """
        filename = filedialog.asksaveasfilename(defaultextension=".session", filetypes=[("Session files", "*.session")])
        if filename:
            session_data = {'files': self.file_entry.get(), 'output': self.output_entry.get(), 'results': self.results}
            save_session(filename, session_data)
            messagebox.showinfo("Success", "Session saved")

    def load_session(self):
        """
        Φορτώνει προηγουμένως αποθηκευμένη συνεδρία από αρχείο.
        Επαναφέρει αρχεία, φάκελο εξόδου και αποτελέσματα.
        """
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
        """
        Δημιουργεί ολοκληρωμένη αναφορά (PDF ή Excel) των αποτελεσμάτων ανάλυσης.
        Ζητά από τον χρήστη τον τύπο αναφοράς που επιθυμεί.
        """
        if not self.results:
            messagebox.showwarning("Warning", "No results to generate report")
            return

        output_dir = self.output_entry.get() or filedialog.askdirectory(title="Select directory for report")
        if not output_dir:
            return

        # Ερώτηση για τύπο αναφοράς
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
        """
        Εμφανίζει τον οδηγό χρήσης της εφαρμογής.
        """
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
        """
        Εμφανίζει πληροφορίες σχετικά με την εφαρμογή.
        """
        about_text = """Professional RNA Analysis Bioinformatics

A comprehensive toolkit for DNA/RNA sequence analysis
designed for research and educational purposes.

Built with Python, BioPython, and tkinter"""
        messagebox.showinfo("About", about_text)

# Σημείο εισόδου της εφαρμογής
if __name__ == "__main__":
    import sys, logging, os
    # Ρύθμιση logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        # Προτίμηση ttkbootstrap Window για αναγνώριση bootstyle
        root = tb.Window(themename="flatly")
        app = EnhancedDNAToolGUI(root)
        logging.info("GUI initialized; entering mainloop")
        root.mainloop()
    except Exception as e:
        logging.exception("Failed to start GUI")
        # Ρητή υπόδειξη για headless περιβάλλον CI/workflow
        if os.name != "nt" and not os.environ.get("DISPLAY"):
            print("Error: No display detected ($DISPLAY unset). GUI cannot launch in this environment.")
        sys.exit(1)
