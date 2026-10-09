import tkinter as tk
from tkinter import messagebox

def nucleotide_profile():
    dna = sequence_entry.get().upper().strip()

    if len(dna) <= 20:
        messagebox.showerror("Invalid sequence", "The DNA sequence must be longer than 20 letters.")
        return

    if any(base not in "ATCG" for base in dna):
        messagebox.showerror("Invalid sequence", "Use only A, T, C and G.")
        return

    total = len(dna)

    result_text.delete("1.0", tk.END)
    result_text.insert(tk.END, f"DNA sequence: {dna}\n")
    result_text.insert(tk.END, f"Length: {total}\n\n")

    for base in "ATCG":
        count = dna.count(base)
        frequency = count / total * 100
        result_text.insert(
            tk.END,
            f"{base} -> {count} occurrences, {frequency:.2f}%\n"
        )

window = tk.Tk()
window.title("DNA Frequency Analyzer")
window.geometry("500x420")
window.resizable(False, False)

title_label = tk.Label(
    window,
    text="DNA Frequency Analyzer",
    font=("Arial", 18, "bold")
)
title_label.pack(pady=20)

instruction_label = tk.Label(
    window,
    text="Enter a DNA sequence longer than 20 letters:"
)
instruction_label.pack()

sequence_entry = tk.Entry(
    window,
    width=50,
    font=("Arial", 12)
)
sequence_entry.pack(pady=10)

analyze_button = tk.Button(
    window,
    text="Analyze Sequence",
    command=nucleotide_profile,
    width=20
)
analyze_button.pack(pady=12)

result_text = tk.Text(
window,
width=50,
height=10,
font=("Consolas", 11)
)
result_text.pack(pady=15)

window.mainloop()
