import psutil

def low_memory_mode():
    return psutil.virtual_memory().total < 4 * 1024 * 1024 * 1024
