import psutil

# Επιστρέφει True αν η συνολική διαθέσιμη RAM είναι κάτω από 4GB
def low_memory_mode():
    # Έλεγχος της συνολικής μνήμης του συστήματος χρησιμοποιώντας τη βιβλιοθήκη psutil
    return psutil.virtual_memory().total < 4 * 1024 * 1024 * 1024  # 4GB σε bytes
