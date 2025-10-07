import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
import re
from Bio.Seq import Seq


def sequence_converter_window(root):
    window = tk.Toplevel(root)
    window.title("Sequence Converter")
    window.geometry("500x400")

    tk.Label(window, text="Input sequence:").pack(pady=5)
    input_text = tk.Text(window, height=5)
    input_text.pack(pady=5, padx=10, fill='x')

    tk.Label(window, text="Output:").pack(pady=5)
    output_text = tk.Text(window, height=10)
    output_text.pack(pady=5, padx=10, fill='both', expand=True)

    def convert():
        sequence = re.sub(r'[^ATGCN]', '', input_text.get(1.0, tk.END).strip().upper())
        if sequence:
            seq_obj = Seq(sequence)
            output = (
                f"Original: {sequence}\n\n"
                f"Reverse: {sequence[::-1]}\n\n"
                f"Complement: {seq_obj.complement()}\n\n"
                f"Reverse Complement: {seq_obj.reverse_complement()}\n\n"
                f"RNA: {seq_obj.transcribe()}\n\n"
            )
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, output)

    tk.Button(window, text="Convert", command=convert).pack(pady=10)


def ask_string(title, prompt, initial=''):
    return simpledialog.askstring(title, prompt, initialvalue=initial)


def ask_integer(title, prompt, initial=100, minval=10, maxval=1000):
    return simpledialog.askinteger(title, prompt, initialvalue=initial, minvalue=minval, maxvalue=maxval)
