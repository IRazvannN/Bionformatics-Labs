import tkinter as tk
from tkinter import messagebox

genetic_code = {
    "UUU": ("Phe", "F"), "UUC": ("Phe", "F"),
    "UUA": ("Leu", "L"), "UUG": ("Leu", "L"),
    "UCU": ("Ser", "S"), "UCC": ("Ser", "S"),
    "UCA": ("Ser", "S"), "UCG": ("Ser", "S"),
    "UAU": ("Tyr", "Y"), "UAC": ("Tyr", "Y"),
    "UAA": ("Stop", "*"), "UAG": ("Stop", "*"),
    "UGU": ("Cys", "C"), "UGC": ("Cys", "C"),
    "UGA": ("Stop", "*"), "UGG": ("Trp", "W"),

    "CUU": ("Leu", "L"), "CUC": ("Leu", "L"),
    "CUA": ("Leu", "L"), "CUG": ("Leu", "L"),
    "CCU": ("Pro", "P"), "CCC": ("Pro", "P"),
    "CCA": ("Pro", "P"), "CCG": ("Pro", "P"),
    "CAU": ("His", "H"), "CAC": ("His", "H"),
    "CAA": ("Gln", "Q"), "CAG": ("Gln", "Q"),
    "CGU": ("Arg", "R"), "CGC": ("Arg", "R"),
    "CGA": ("Arg", "R"), "CGG": ("Arg", "R"),

    "AUU": ("Ile", "I"), "AUC": ("Ile", "I"),
    "AUA": ("Ile", "I"), "AUG": ("Met", "M"),
    "ACU": ("Thr", "T"), "ACC": ("Thr", "T"),
    "ACA": ("Thr", "T"), "ACG": ("Thr", "T"),
    "AAU": ("Asn", "N"), "AAC": ("Asn", "N"),
    "AAA": ("Lys", "K"), "AAG": ("Lys", "K"),
    "AGU": ("Ser", "S"), "AGC": ("Ser", "S"),
    "AGA": ("Arg", "R"), "AGG": ("Arg", "R"),

    "GUU": ("Val", "V"), "GUC": ("Val", "V"),
    "GUA": ("Val", "V"), "GUG": ("Val", "V"),
    "GCU": ("Ala", "A"), "GCC": ("Ala", "A"),
    "GCA": ("Ala", "A"), "GCG": ("Ala", "A"),
    "GAU": ("Asp", "D"), "GAC": ("Asp", "D"),
    "GAA": ("Glu", "E"), "GAG": ("Glu", "E"),
    "GGU": ("Gly", "G"), "GGC": ("Gly", "G"),
    "GGA": ("Gly", "G"), "GGG": ("Gly", "G")
}

def translate_sequence():
    dna = sequence_entry.get("1.0", tk.END).upper().replace("\n", "").replace(" ", "")

    if not dna:
        messagebox.showerror(
            "Invalid sequence",
            "Enter a DNA sequence."
        )
        return

    if any(base not in "ATCG" for base in dna):
        messagebox.showerror(
            "Invalid sequence",
            "Use only A, T, C and G."
        )
        return

    if len(dna) < 100:
        messagebox.showwarning(
            "Short sequence",
            "The sequence contains fewer than 100 letters."
        )

    rna = dna.replace("T", "U")

    start = rna.find("AUG")

    if start == -1:
        messagebox.showerror(
            "No start codon",
            "No AUG start codon was found in the sequence."
        )
        return

    amino_three = []
    amino_one = []
    codons = []
    stop_codon = None

    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        codons.append(codon)

        amino = genetic_code[codon]

        if amino[0] == "Stop":
            stop_codon = codon
            break

        amino_three.append(amino[0])
        amino_one.append(amino[1])

    result_text.delete("1.0", tk.END)

    result_text.insert(
        tk.END,
        f"DNA length: {len(dna)}\n\n"
    )

    result_text.insert(
        tk.END,
        "Transcription - RNA sequence:\n"
    )
    result_text.insert(
        tk.END,
        rna + "\n\n"
    )

    result_text.insert(
        tk.END,
        f"First AUG found at position: {start + 1}\n\n"
    )

    result_text.insert(
        tk.END,
        "Reading frame:\n"
    )
    result_text.insert(
        tk.END,
        "-".join(codons) + "\n\n"
    )

    result_text.insert(
        tk.END,
        "Translation - 3-letter amino acid sequence:\n"
    )
    result_text.insert(
        tk.END,
        "-".join(amino_three) + "\n\n"
    )

    result_text.insert(
        tk.END,
        "Translation - 1-letter amino acid sequence:\n"
    )
    result_text.insert(
        tk.END,
        "-".join(amino_one) + "\n\n"
    )

    if stop_codon:
        result_text.insert(
            tk.END,
            f"Translation stopped at: {stop_codon}"
        )
    else:
        result_text.insert(
            tk.END,
            "No stop codon was found in this reading frame."
        )


window = tk.Tk()
window.title("DNA Transcription and Translation")
window.geometry("760x700")
window.resizable(False, False)

title_label = tk.Label(
    window,
    text="DNA Transcription and Translation",
    font=("Arial", 18, "bold")
)
title_label.pack(pady=20)

instruction_label = tk.Label(
    window,
    text="Enter a DNA sequence:"
)
instruction_label.pack()

sequence_entry = tk.Text(
    window,
    width=80,
    height=6,
    font=("Consolas", 11)
)
sequence_entry.pack(pady=10)

translate_button = tk.Button(
    window,
    text="Translate Sequence",
    command=translate_sequence,
    width=22
)
translate_button.pack(pady=10)

result_text = tk.Text(
    window,
    width=85,
    height=28,
    font=("Consolas", 10)
)
result_text.pack(pady=15)

window.mainloop()