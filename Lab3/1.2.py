import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk
from itertools import product
from urllib.request import urlopen, Request
from urllib.parse import quote
import zipfile
import io
import os

codon_counts = {}


def create_codon_dictionary():
    bases = "ATCG"
    counts = {}

    for letters in product(bases, repeat=3):
        codon = "".join(letters)
        counts[codon] = 0

    return counts


def fasta_to_sequence(fasta):
    sequence = ""

    for line in fasta.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            continue

        sequence += line.upper()

    return sequence


def read_fasta(file_path):
    with open(file_path, "r") as file:
        fasta = file.read()

    return fasta_to_sequence(fasta)


def count_codons(sequence):
    counts = create_codon_dictionary()

    for i in range(len(sequence) - 2):
        codon = sequence[i:i + 3]

        if codon in counts:
            counts[codon] += 1

    return counts


def analyze_sequence(sequence, name):
    global codon_counts

    if not sequence:
        messagebox.showerror(
            "Error",
            "No sequence was found."
        )
        return

    valid_bases = sum(
        1 for base in sequence
        if base in "ATCG"
    )

    if valid_bases == 0:
        messagebox.showerror(
            "Error",
            "Invalid DNA sequence."
        )
        return

    codon_counts = count_codons(sequence)

    file_label.config(
        text=name
    )

    length_label.config(
        text=f"Sequence length: {len(sequence)}"
    )

    show_top_10()
    show_all_codons()


def open_fasta():
    file_path = filedialog.askopenfilename(
        title="Open File",
        filetypes=[
            ("FASTA files", "*.fna *.fasta *.fa"),
            ("All files", "*.*")
        ]
    )

    if not file_path:
        return

    try:
        sequence = read_fasta(file_path)

        analyze_sequence(
            sequence,
            os.path.basename(file_path)
        )

        status_label.config(
            text="File loaded"
        )

    except Exception as error:
        messagebox.showerror(
            "Error",
            str(error)
        )


def get_ncbi_sequence():
    accession = accession_entry.get().strip().upper()

    if not accession:
        messagebox.showerror(
            "Error",
            "Enter a GCF accession."
        )
        return

    if not accession.startswith("GCF_"):
        messagebox.showerror(
            "Error",
            "Enter an accession such as GCF_000859625.1"
        )
        return

    try:
        status_label.config(
            text="Loading..."
        )

        window.update()

        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + quote(accession)
            + "/download?include_annotation_type=GENOME_FASTA"
        )

        request = Request(
            url,
            headers={
                "User-Agent": "Python Codon Frequency Application"
            }
        )

        with urlopen(request, timeout=30) as response:
            data = response.read()

        zip_file = zipfile.ZipFile(
            io.BytesIO(data)
        )

        fasta_files = [
            name for name in zip_file.namelist()
            if name.endswith(".fna")
        ]

        if not fasta_files:
            messagebox.showerror(
                "Error",
                "No genomic FASTA was found."
            )

            status_label.config(
                text=""
            )

            return

        fasta = ""

        for fasta_file in fasta_files:
            fasta += zip_file.read(
                fasta_file
            ).decode("utf-8")

        sequence = fasta_to_sequence(fasta)

        analyze_sequence(
            sequence,
            accession
        )

        status_label.config(
            text="Sequence loaded"
        )

    except Exception as error:
        status_label.config(
            text=""
        )

        messagebox.showerror(
            "Error",
            f"Could not retrieve the sequence.\n\n{error}"
        )


def show_top_10():
    top_text.delete(
        "1.0",
        tk.END
    )

    sorted_codons = sorted(
        codon_counts.items(),
        key=lambda item: item[1],
        reverse=True
    )

    top_10 = sorted_codons[:10]

    total = sum(
        codon_counts.values()
    )

    top_text.insert(
        tk.END,
        "Top 10\n\n"
    )

    for position, (codon, count) in enumerate(
        top_10,
        start=1
    ):
        if total > 0:
            frequency = (
                count / total
            ) * 100
        else:
            frequency = 0

        top_text.insert(
            tk.END,
            f"{position}.  {codon}    "
            f"Count: {count}    "
            f"Frequency: {frequency:.2f}%\n"
        )


def show_all_codons():
    for item in codon_table.get_children():
        codon_table.delete(item)

    sorted_codons = sorted(
        codon_counts.items(),
        key=lambda item: item[1],
        reverse=True
    )

    total = sum(
        codon_counts.values()
    )

    for position, (codon, count) in enumerate(
        sorted_codons,
        start=1
    ):
        if total > 0:
            frequency = (
                count / total
            ) * 100
        else:
            frequency = 0

        codon_table.insert(
            "",
            "end",
            values=(
                position,
                codon,
                count,
                f"{frequency:.2f}%"
            )
        )


window = tk.Tk()

window.title("Codon Frequency")

window.geometry("900x750")
window.minsize(750, 600)

title_label = tk.Label(
    window,
    text="Codon Frequency",
    font=("Arial", 20, "bold")
)

title_label.pack(
    pady=20
)

ncbi_frame = tk.Frame(
    window
)

ncbi_frame.pack(
    pady=5
)

accession_label = tk.Label(
    ncbi_frame,
    text="NCBI accession:"
)

accession_label.pack(
    side="left",
    padx=5
)

accession_entry = tk.Entry(
    ncbi_frame,
    width=25
)

accession_entry.pack(
    side="left",
    padx=5
)

get_button = tk.Button(
    ncbi_frame,
    text="Get Sequence",
    command=get_ncbi_sequence,
    width=15
)

get_button.pack(
    side="left",
    padx=5
)

status_label = tk.Label(
    window,
    text=""
)

status_label.pack(
    pady=3
)

or_label = tk.Label(
    window,
    text="or"
)

or_label.pack(
    pady=3
)

open_button = tk.Button(
    window,
    text="Open File",
    command=open_fasta,
    width=20
)

open_button.pack(
    pady=7
)

file_label = tk.Label(
    window,
    text="No sequence selected"
)

file_label.pack(
    pady=5
)

length_label = tk.Label(
    window,
    text="Sequence length: -"
)

length_label.pack(
    pady=3
)

notebook = ttk.Notebook(
    window
)

notebook.pack(
    fill="both",
    expand=True,
    padx=15,
    pady=15
)

top_tab = tk.Frame(
    notebook
)

all_tab = tk.Frame(
    notebook
)

notebook.add(
    top_tab,
    text="Top 10"
)

notebook.add(
    all_tab,
    text="All Codons"
)

top_text = tk.Text(
    top_tab,
    font=("Consolas", 12),
    wrap="word"
)

top_text.pack(
    fill="both",
    expand=True,
    padx=10,
    pady=10
)

table_frame = tk.Frame(
    all_tab
)

table_frame.pack(
    fill="both",
    expand=True,
    padx=10,
    pady=10
)

scrollbar = ttk.Scrollbar(
    table_frame,
    orient="vertical"
)

scrollbar.pack(
    side="right",
    fill="y"
)

codon_table = ttk.Treeview(
    table_frame,
    columns=(
        "rank",
        "codon",
        "count",
        "frequency"
    ),
    show="headings",
    yscrollcommand=scrollbar.set
)

codon_table.heading(
    "rank",
    text="Rank"
)

codon_table.heading(
    "codon",
    text="Codon"
)

codon_table.heading(
    "count",
    text="Count"
)

codon_table.heading(
    "frequency",
    text="Frequency"
)

codon_table.column(
    "rank",
    width=100,
    anchor="center"
)

codon_table.column(
    "codon",
    width=150,
    anchor="center"
)

codon_table.column(
    "count",
    width=150,
    anchor="center"
)

codon_table.column(
    "frequency",
    width=150,
    anchor="center"
)

codon_table.pack(
    side="left",
    fill="both",
    expand=True
)

scrollbar.config(
    command=codon_table.yview
)

window.mainloop()