import tkinter as tk
from tkinter import messagebox

def combination_profile():
    dna = sequence_entry.get().upper().strip()

    if len(dna) <= 20:
        messagebox.showerror(
            "Invalid sequence",
            "The DNA sequence must be longer than 20 letters."
        )
        return

    if any(base not in "ATCG" for base in dna):
        messagebox.showerror(
            "Invalid sequence",
            "Use only A, T, C and G."
        )
        return

    pairs = {}
    triplets = {}

    for i in range(len(dna) - 1):
        combination = dna[i:i + 2]

        if combination in pairs:
            pairs[combination] += 1
        else:
            pairs[combination] = 1

    for i in range(len(dna) - 2):
        combination = dna[i:i + 3]

        if combination in triplets:
            triplets[combination] += 1
        else:
            triplets[combination] = 1

    total_pairs = len(dna) - 1
    total_triplets = len(dna) - 2

    result_text.delete("1.0", tk.END)

    result_text.insert(
        tk.END,
        f"DNA sequence: {dna}\n"
    )

    result_text.insert(
        tk.END,
        f"Length: {len(dna)}\n\n"
    )

    result_text.insert(
        tk.END,
        "2-letter combinations:\n"
    )

    for combination, count in sorted(pairs.items()):
        frequency = count / total_pairs * 100

        result_text.insert(
            tk.END,
            f"{combination} -> {count} occurrences, {frequency:.2f}%\n"
        )

    result_text.insert(
        tk.END,
        "\n3-letter combinations:\n"
    )

    for combination, count in sorted(triplets.items()):
        frequency = count / total_triplets * 100

        result_text.insert(
            tk.END,
            f"{combination} -> {count} occurrences, {frequency:.2f}%\n"
        )


window = tk.Tk()
window.title("DNA Combination Analyzer")
window.geometry("620x650")
window.resizable(False, False)

title_label = tk.Label(
    window,
    text="DNA Combination Analyzer",
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
    width=55,
    font=("Arial", 12)
)
sequence_entry.pack(pady=10)

analyze_button = tk.Button(
    window,
    text="Analyze Combinations",
    command=combination_profile,
    width=22
)
analyze_button.pack(pady=12)

result_text = tk.Text(
    window,
    width=65,
    height=28,
    font=("Consolas", 10)
)
result_text.pack(pady=10)

window.mainloop()