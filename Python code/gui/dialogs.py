
"""
Μονάδα: Dialogs
Παρέχει διάφορα παράθυρα διαλόγου για την εφαρμογή βιοπληροφορικής.
Περιλαμβάνει μετατροπέα αλληλουχιών και βοηθητικές συναρτήσεις για είσοδο δεδομένων.
"""
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
import re
from Bio.Seq import Seq


def sequence_converter_window(root):
    """
    Δημιουργεί παράθυρο διαλόγου για μετατροπή αλληλουχιών DNA.
    Επιτρέπει στον χρήστη να εισάγει μια αλληλουχία και να δει:
    - Την αντίστροφη αλληλουχία
    - Το συμπλήρωμα
    - Το αντίστροφο συμπλήρωμα
    - Τη μεταγραφή σε RNA
    
    Παράμετροι:
        root (Tk): Το κύριο παράθυρο της εφαρμογής
    """
    # Δημιουργία νέου παραθύρου
    window = tk.Toplevel(root)
    window.title("Sequence Converter")
    window.geometry("500x400")

    # Ετικέτα για το πεδίο εισαγωγής
    tk.Label(window, text="Input sequence:").pack(pady=5)
    
    # Πεδίο κειμένου για εισαγωγή αλληλουχίας
    input_text = tk.Text(window, height=5)
    input_text.pack(pady=5, padx=10, fill='x')

    # Ετικέτα για το πεδίο εξόδου
    tk.Label(window, text="Output:").pack(pady=5)
    
    # Πεδίο κειμένου για εμφάνιση αποτελεσμάτων
    output_text = tk.Text(window, height=10)
    output_text.pack(pady=5, padx=10, fill='both', expand=True)

    def convert():
        """
        Μετατρέπει την εισαγόμενη αλληλουχία DNA σε διάφορες μορφές.
        Καθαρίζει την είσοδο από μη έγκυρους χαρακτήρες και εμφανίζει:
        - Αρχική αλληλουχία
        - Αντίστροφη
        - Συμπλήρωμα
        - Αντίστροφο συμπλήρωμα
        - Μεταγραφή σε RNA
        """
        # Καθαρισμός αλληλουχίας: κρατάμε μόνο A, T, G, C, N
        sequence = re.sub(r'[^ATGCN]', '', input_text.get(1.0, tk.END).strip().upper())
        
        if sequence:
            # Δημιουργία αντικειμένου BioPython Seq
            seq_obj = Seq(sequence)
            
            # Δημιουργία κειμένου εξόδου με όλες τις μετατροπές
            output = (
                f"Original: {sequence}\n\n"
                f"Reverse: {sequence[::-1]}\n\n"
                f"Complement: {seq_obj.complement()}\n\n"
                f"Reverse Complement: {seq_obj.reverse_complement()}\n\n"
                f"RNA: {seq_obj.transcribe()}\n\n"
            )
            
            # Καθαρισμός και εμφάνιση αποτελεσμάτων
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, output)

    # Κουμπί για εκτέλεση μετατροπής
    tk.Button(window, text="Convert", command=convert).pack(pady=10)


def ask_string(title, prompt, initial=''):
    """
    Εμφανίζει διάλογο για εισαγωγή κειμένου από τον χρήστη.
    
    Παράμετροι:
        title (str): Τίτλος του παραθύρου διαλόγου
        prompt (str): Μήνυμα προτροπής για τον χρήστη
        initial (str): Αρχική τιμή (προεπιλογή: κενό string)
        
    Επιστρέφει:
        str: Το κείμενο που εισήγαγε ο χρήστης ή None αν ακυρώθηκε
    """
    return simpledialog.askstring(title, prompt, initialvalue=initial)


def ask_integer(title, prompt, initial=100, minval=10, maxval=1000):
    """
    Εμφανίζει διάλογο για εισαγωγή ακέραιου αριθμού από τον χρήστη.
    
    Παράμετροι:
        title (str): Τίτλος του παραθύρου διαλόγου
        prompt (str): Μήνυμα προτροπής για τον χρήστη
        initial (int): Αρχική τιμή (προεπιλογή: 100)
        minval (int): Ελάχιστη επιτρεπτή τιμή (προεπιλογή: 10)
        maxval (int): Μέγιστη επιτρεπτή τιμή (προεπιλογή: 1000)
        
    Επιστρέφει:
        int: Ο ακέραιος που εισήγαγε ο χρήστης ή None αν ακυρώθηκε
    """
    return simpledialog.askinteger(title, prompt, initialvalue=initial, minvalue=minval, maxvalue=maxval)
