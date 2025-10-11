
"""
Μονάδα: Εφαρμογή
Κύριο σημείο εισόδου για την εφαρμογή βιοπληροφορικής DNA/RNA.
Αρχικοποιεί το GUI και εκκινεί τον κύριο βρόχο της εφαρμογής.
"""

from tabs import EnhancedDNAToolGUI
import ttkbootstrap as tb

def main():
        """
        Κύρια συνάρτηση εκκίνησης της εφαρμογής.
        Δημιουργεί το κύριο παράθυρο με θέμα 'superhero' και 
        αρχικοποιεί το GUI της εφαρμογής.
        """
        # Δημιουργία κύριου παραθύρου με θέμα ttkbootstrap
        root = tb.Window(themename='superhero')

        # Αρχικοποίηση του GUI της εφαρμογής
        app = EnhancedDNAToolGUI(root)

        # Εκκίνηση του κύριου βρόχου της εφαρμογής
        root.mainloop()

if __name__ == '__main__':
        # Εκτέλεση της κύριας συνάρτησης όταν το script εκτελείται απευθείας
        main()
