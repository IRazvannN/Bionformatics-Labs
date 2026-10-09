import tkinter as tk
from tkinter import messagebox
from tkinter import ttk

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

stop_codons = {"UAA", "UAG", "UGA"}

def clean_sequence():
    raw = sequence_entry.get("1.0", tk.END).strip()
    lines = raw.splitlines()

    sequence_lines = []

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            continue

        sequence_lines.append(line)

    dna = "".join(sequence_lines).upper().replace(" ", "")

    return dna

def analyze_sequence():
    dna = clean_sequence()

    if not dna:
        messagebox.showerror(
            "Invalid sequence",
            "Enter a DNA sequence."
        )
        return

    if any(base not in "ATCG" for base in dna):
        messagebox.showerror(
            "Invalid sequence",
            "The sequence must contain only A, T, C and G."
        )
        return

    rna = dna.replace("T", "U")

    products = []
    aug_count = 0

    for start in range(len(rna) - 2):
        if rna[start:start + 3] != "AUG":
            continue

        aug_count += 1

        codons = []
        amino_three = []
        amino_one = []
        stop_codon = None
        stop_position = None

        for i in range(start, len(rna) - 2, 3):
            codon = rna[i:i + 3]
            codons.append(codon)

            if codon in stop_codons:
                stop_codon = codon
                stop_position = i
                break

            amino = genetic_code[codon]
            amino_three.append(amino[0])
            amino_one.append(amino[1])

        if stop_codon:
            products.append({
                "start": start,
                "stop": stop_position,
                "stop_codon": stop_codon,
                "codons": codons,
                "three": amino_three,
                "one": amino_one
            })

    result_text.delete("1.0", tk.END)

    result_text.insert(
        tk.END,
        f"DNA length: {len(dna)}\n"
    )

    result_text.insert(
        tk.END,
        f"Total AUG codons found: {aug_count}\n"
    )

    result_text.insert(
        tk.END,
        f"Possible translation products: {len(products)}\n\n"
    )

    if not products:
        result_text.insert(
            tk.END,
            "No valid AUG-to-stop translation products were found."
        )
        return

    for number, product in enumerate(products, start=1):
        result_text.insert(
            tk.END,
            f"Product {number}\n"
        )

        result_text.insert(
            tk.END,
            f"Start AUG position: {product['start'] + 1}\n"
        )

        result_text.insert(
            tk.END,
            f"Stop codon: {product['stop_codon']} at position {product['stop'] + 1}\n"
        )

        result_text.insert(
            tk.END,
            f"Protein length: {len(product['one'])} amino acids\n"
        )

        result_text.insert(
            tk.END,
            "Codons:\n"
        )

        result_text.insert(
            tk.END,
            "-".join(product["codons"]) + "\n"
        )

        result_text.insert(
            tk.END,
            "3-letter sequence:\n"
        )

        result_text.insert(
            tk.END,
            "-".join(product["three"]) + "\n"
        )

        result_text.insert(
            tk.END,
            "1-letter sequence:\n"
        )

        result_text.insert(
            tk.END,
            "-".join(product["one"]) + "\n"
        )

        result_text.insert(
            tk.END,
            "\n" + "=" * 70 + "\n\n"
        )

def toggle_fullscreen(event=None):
    window.attributes("-fullscreen", not window.attributes("-fullscreen"))

def exit_fullscreen(event=None):
    window.attributes("-fullscreen", False)

window = tk.Tk()
window.title("Influenza Translation Product Analyzer")
window.geometry("1100x800")
window.minsize(800, 600)
window.resizable(True, True)

window.bind("<F11>", toggle_fullscreen)
window.bind("<Escape>", exit_fullscreen)

main_frame = tk.Frame(window)
main_frame.pack(fill="both", expand=True, padx=20, pady=15)

title_label = tk.Label(
    main_frame,
    text="Influenza Translation Product Analyzer",
    font=("Arial", 18, "bold")
)
title_label.pack(pady=10)

instruction_label = tk.Label(
    main_frame,
    text="Paste a DNA sequence or FASTA record:"
)
instruction_label.pack()

input_frame = tk.Frame(main_frame)
input_frame.pack(fill="x", pady=10)

input_scrollbar = ttk.Scrollbar(
    input_frame,
    orient="vertical"
)
input_scrollbar.pack(side="right", fill="y")

sequence_entry = tk.Text(
    input_frame,
    height=10,
    font=("Consolas", 10),
    wrap="word",
    yscrollcommand=input_scrollbar.set
)
sequence_entry.pack(
    side="left",
    fill="both",
    expand=True
)

input_scrollbar.config(
    command=sequence_entry.yview
)

analyze_button = tk.Button(
    main_frame,
    text="Find Translation Products",
    command=analyze_sequence,
    width=25
)
analyze_button.pack(pady=10)

output_frame = tk.Frame(main_frame)
output_frame.pack(
    fill="both",
    expand=True,
    pady=10
)

output_scrollbar = ttk.Scrollbar(
    output_frame,
    orient="vertical"
)
output_scrollbar.pack(
    side="right",
    fill="y"
)

result_text = tk.Text(
    output_frame,
    font=("Consolas", 10),
    wrap="word",
    yscrollcommand=output_scrollbar.set
)
result_text.pack(
    side="left",
    fill="both",
    expand=True
)

output_scrollbar.config(
    command=result_text.yview
)

window.mainloop()