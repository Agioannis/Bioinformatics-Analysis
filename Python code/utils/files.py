
"""
Μονάδα: Διαχείριση Αρχείων
Παρέχει λειτουργίες για αποθήκευση και φόρτωση συνεδριών ανάλυσης.
Χρησιμοποιεί το module pickle για σειριοποίηση δεδομένων.
"""

import pickle
import logging
from typing import Dict, Optional

def save_session(filename: str, data: Dict) -> bool:
    """
    Αποθηκεύει μια συνεδρία ανάλυσης σε αρχείο.
    
    Παράμετροι:
        filename (str): Το όνομα του αρχείου για αποθήκευση
        data (Dict): Το λεξικό με τα δεδομένα της συνεδρίας
    
    Επιστρέφει:
        bool: True αν η αποθήκευση ήταν επιτυχής, False αλλιώς
    """
    try:
        # Άνοιγμα αρχείου σε binary mode για εγγραφή
        with open(filename, 'wb') as f:
            # Σειριοποίηση και αποθήκευση δεδομένων με pickle
            pickle.dump(data, f)
        return True
    except Exception as e:
        # Καταγραφή σφάλματος στο log
        logging.error(f"Save session failed: {e}")
        return False

def load_session(filename: str) -> Optional[Dict]:
    """
    Φορτώνει μια αποθηκευμένη συνεδρία ανάλυσης από αρχείο.
    
    Παράμετροι:
        filename (str): Το όνομα του αρχείου προς φόρτωση
    
    Επιστρέφει:
        Optional[Dict]: Το λεξικό με τα δεδομένα της συνεδρίας αν η φόρτωση 
                       ήταν επιτυχής, None αλλιώς
    """
    try:
        # Άνοιγμα αρχείου σε binary mode για ανάγνωση
        with open(filename, 'rb') as f:
            # Αποσειριοποίηση και επιστροφή δεδομένων
            return pickle.load(f)
    except Exception as e:
        # Καταγραφή σφάλματος στο log
        logging.error(f"Load session failed: {e}")
        return None
