import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk
from itertools import product
from urllib.request import urlopen, Request
from urllib.parse import quote, urlencode
from zipfile import ZipFile
from io import BytesIO
import json
import os

codon_counts = {}
USDA_API_KEY = "DEMO_KEY"

genetic_code = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine",
    "TTA": "Leucine", "TTG": "Leucine",
    "TCT": "Serine", "TCC": "Serine",
    "TCA": "Serine", "TCG": "Serine",
    "TAT": "Tyrosine", "TAC": "Tyrosine",
    "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cysteine", "TGC": "Cysteine",
    "TGA": "Stop", "TGG": "Tryptophan",

    "CTT": "Leucine", "CTC": "Leucine",
    "CTA": "Leucine", "CTG": "Leucine",
    "CCT": "Proline", "CCC": "Proline",
    "CCA": "Proline", "CCG": "Proline",
    "CAT": "Histidine", "CAC": "Histidine",
    "CAA": "Glutamine", "CAG": "Glutamine",
    "CGT": "Arginine", "CGC": "Arginine",
    "CGA": "Arginine", "CGG": "Arginine",

    "ATT": "Isoleucine", "ATC": "Isoleucine",
    "ATA": "Isoleucine", "ATG": "Methionine",
    "ACT": "Threonine", "ACC": "Threonine",
    "ACA": "Threonine", "ACG": "Threonine",
    "AAT": "Asparagine", "AAC": "Asparagine",
    "AAA": "Lysine", "AAG": "Lysine",
    "AGT": "Serine", "AGC": "Serine",
    "AGA": "Arginine", "AGG": "Arginine",

    "GTT": "Valine", "GTC": "Valine",
    "GTA": "Valine", "GTG": "Valine",
    "GCT": "Alanine", "GCC": "Alanine",
    "GCA": "Alanine", "GCG": "Alanine",
    "GAT": "Aspartic acid", "GAC": "Aspartic acid",
    "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGT": "Glycine", "GGC": "Glycine",
    "GGA": "Glycine", "GGG": "Glycine"
}


def create_codon_dictionary():
    counts = {}

    for letters in product("ATCG", repeat=3):
        counts["".join(letters)] = 0

    return counts


def fasta_to_sequence(fasta):
    sequence = ""

    for line in fasta.splitlines():
        line = line.strip()

        if line and not line.startswith(">"):
            sequence += line.upper()

    return sequence


def read_fasta(file_path):
    with open(file_path, "r") as file:
        return fasta_to_sequence(file.read())


def count_codons(sequence):
    counts = create_codon_dictionary()

    for i in range(len(sequence) - 2):
        codon = sequence[i:i + 3]

        if codon in counts:
            counts[codon] += 1

    return counts


def get_sorted_codons():
    return sorted(
        codon_counts.items(),
        key=lambda item: item[1],
        reverse=True
    )


def get_top_amino_acids():
    amino_acids = []

    for codon, count in get_sorted_codons()[:10]:
        amino = genetic_code[codon]

        if amino != "Stop" and amino not in amino_acids:
            amino_acids.append(amino)

    return amino_acids


def analyze_sequence(sequence, name):
    global codon_counts

    if not sequence:
        messagebox.showerror("Error", "No sequence was found.")
        return

    valid = sum(1 for base in sequence if base in "ATCG")

    if valid == 0:
        messagebox.showerror("Error", "Invalid DNA sequence.")
        return

    codon_counts = count_codons(sequence)

    file_label.config(text=name)
    length_label.config(text=f"Sequence length: {len(sequence)}")

    show_top_10()
    show_all_codons()
    prepare_diet()


def get_ncbi_sequence():
    accession = accession_entry.get().strip().upper()

    if not accession:
        messagebox.showerror(
            "Error",
            "Enter an assembly accession."
        )
        return

    if not accession.startswith(("GCF_", "GCA_")):
        messagebox.showerror(
            "Error",
            "Enter a GCF_ or GCA_ assembly accession."
        )
        return

    try:
        status_label.config(text="Loading genome...")
        window.update()

        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/"
            "genome/accession/"
            + quote(accession)
            + "/download?include_annotation_type=GENOME_FASTA"
        )

        request = Request(
            url,
            headers={
                "User-Agent": "BioinformaticsStudentApp/1.0"
            }
        )

        with urlopen(request, timeout=30) as response:
            zip_data = response.read()

        with ZipFile(BytesIO(zip_data)) as archive:
            fasta_files = [
                name for name in archive.namelist()
                if name.endswith("_genomic.fna")
            ]

            if not fasta_files:
                raise Exception("No genomic FASTA was found.")

            fasta = archive.read(
                fasta_files[0]
            ).decode("utf-8")

        sequence = fasta_to_sequence(fasta)

        analyze_sequence(sequence, accession)

        status_label.config(text="Sequence loaded")

    except Exception as error:
        status_label.config(text="")

        messagebox.showerror(
            "Error",
            "Could not retrieve the sequence.\n\n"
            + str(error)
        )


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

    except Exception as error:
        messagebox.showerror("Error", str(error))


def show_top_10():
    for item in top_table.get_children():
        top_table.delete(item)

    top_10 = get_sorted_codons()[:10]
    total = sum(codon_counts.values())

    for position, (codon, count) in enumerate(top_10, start=1):
        frequency = count / total * 100 if total else 0

        top_table.insert(
            "",
            "end",
            values=(
                position,
                codon,
                count,
                f"{frequency:.2f}%",
                genetic_code[codon]
            )
        )


def show_all_codons():
    for item in codon_table.get_children():
        codon_table.delete(item)

    total = sum(codon_counts.values())

    for position, (codon, count) in enumerate(
        get_sorted_codons(),
        start=1
    ):
        frequency = count / total * 100 if total else 0

        codon_table.insert(
            "",
            "end",
            values=(
                position,
                codon,
                count,
                f"{frequency:.2f}%",
                genetic_code[codon]
            )
        )


def prepare_diet():
    diet_text.delete("1.0", tk.END)

    amino_acids = get_top_amino_acids()

    diet_text.insert(
        tk.END,
        "Amino acids from Top 10 codons:\n\n"
    )

    for amino in amino_acids:
        diet_text.insert(
            tk.END,
            f"- {amino}\n"
        )

    diet_text.insert(
        tk.END,
        "\nPress 'Find Foods' to search USDA FoodData Central.\n"
    )


def search_usda_foods():
    if not codon_counts:
        messagebox.showerror(
            "Error",
            "Load a genome first."
        )
        return

    amino_acids = get_top_amino_acids()

    try:
        diet_status.config(text="Searching USDA...")
        window.update()

        query = urlencode({
            "api_key": USDA_API_KEY,
            "query": "food",
            "pageSize": 100
        })

        url = (
            "https://api.nal.usda.gov/fdc/v1/foods/search?"
            + query
        )

        request = Request(
            url,
            headers={
                "User-Agent": "BioinformaticsStudentApp/1.0"
            }
        )

        with urlopen(request, timeout=30) as response:
            data = json.loads(
                response.read().decode("utf-8")
            )

        foods = data.get("foods", [])

        results = []

        for food in foods:
            nutrients = food.get("foodNutrients", [])

            amino_values = {}
            protein = None

            for nutrient in nutrients:
                name = nutrient.get("nutrientName", "")
                value = nutrient.get("value")

                if value is None:
                    continue

                if name.lower() == "protein":
                    protein = value

                for amino in amino_acids:
                    if amino.lower() in name.lower():
                        amino_values[amino] = value

            if len(amino_values) == len(amino_acids):
                total_amino = sum(amino_values.values())

                results.append({
                    "name": food.get(
                        "description",
                        "Unknown food"
                    ),
                    "fdc_id": food.get("fdcId"),
                    "protein": protein,
                    "total": total_amino,
                    "values": amino_values
                })

        results.sort(
            key=lambda food: food["total"]
        )

        show_food_results(
            results[:10],
            amino_acids
        )

        diet_status.config(
            text="Source: USDA FoodData Central"
        )

    except Exception as error:
        diet_status.config(text="")

        messagebox.showerror(
            "USDA Error",
            "Could not retrieve food information.\n\n"
            + str(error)
        )


def show_food_results(results, amino_acids):
    diet_text.delete("1.0", tk.END)

    diet_text.insert(
        tk.END,
        "Amino acids from Top 10 codons:\n\n"
    )

    for amino in amino_acids:
        diet_text.insert(
            tk.END,
            f"- {amino}\n"
        )

    diet_text.insert(
        tk.END,
        "\nLowest matching foods found in USDA data:\n\n"
    )

    if not results:
        diet_text.insert(
            tk.END,
            "No foods containing measurements for all "
            "selected amino acids were found.\n"
        )
        return

    for position, food in enumerate(results, start=1):
        diet_text.insert(
            tk.END,
            f"{position}. {food['name']}\n"
        )

        diet_text.insert(
            tk.END,
            f"   FDC ID: {food['fdc_id']}\n"
        )

        if food["protein"] is not None:
            diet_text.insert(
                tk.END,
                f"   Protein: {food['protein']} g\n"
            )

        diet_text.insert(
            tk.END,
            "   Amino acids:\n"
        )

        for amino in amino_acids:
            value = food["values"][amino]

            diet_text.insert(
                tk.END,
                f"      {amino}: {value} g\n"
            )

        diet_text.insert(
            tk.END,
            f"   Combined: {food['total']:.3f} g\n\n"
        )

    diet_text.insert(
        tk.END,
        "Foods are ranked using amino-acid values "
        "reported by USDA FoodData Central.\n"
    )


window = tk.Tk()

window.title("Codon Frequency")
window.geometry("1050x750")
window.minsize(850, 600)

title_label = tk.Label(
    window,
    text="Codon Frequency",
    font=("Arial", 20, "bold")
)

title_label.pack(pady=20)

ncbi_frame = tk.Frame(window)
ncbi_frame.pack(pady=5)

accession_label = tk.Label(
    ncbi_frame,
    text="Assembly:"
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

status_label.pack(pady=3)

or_label = tk.Label(
    window,
    text="or"
)

or_label.pack(pady=3)

open_button = tk.Button(
    window,
    text="Open File",
    command=open_fasta,
    width=20
)

open_button.pack(pady=7)

file_label = tk.Label(
    window,
    text="No sequence selected"
)

file_label.pack(pady=5)

length_label = tk.Label(
    window,
    text="Sequence length: -"
)

length_label.pack(pady=3)

notebook = ttk.Notebook(window)

notebook.pack(
    fill="both",
    expand=True,
    padx=15,
    pady=15
)

top_tab = tk.Frame(notebook)
all_tab = tk.Frame(notebook)
diet_tab = tk.Frame(notebook)

notebook.add(top_tab, text="Top 10")
notebook.add(all_tab, text="All Codons")
notebook.add(diet_tab, text="Diet")

top_table = ttk.Treeview(
    top_tab,
    columns=(
        "rank",
        "codon",
        "count",
        "frequency",
        "amino"
    ),
    show="headings"
)

top_table.heading("rank", text="Rank")
top_table.heading("codon", text="Codon")
top_table.heading("count", text="Count")
top_table.heading("frequency", text="Frequency")
top_table.heading("amino", text="Amino Acid")

top_table.column("rank", width=80, anchor="center")
top_table.column("codon", width=100, anchor="center")
top_table.column("count", width=100, anchor="center")
top_table.column("frequency", width=120, anchor="center")
top_table.column("amino", width=180, anchor="center")

top_table.pack(
    fill="both",
    expand=True,
    padx=10,
    pady=10
)

table_frame = tk.Frame(all_tab)

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
        "frequency",
        "amino"
    ),
    show="headings",
    yscrollcommand=scrollbar.set
)

codon_table.heading("rank", text="Rank")
codon_table.heading("codon", text="Codon")
codon_table.heading("count", text="Count")
codon_table.heading("frequency", text="Frequency")
codon_table.heading("amino", text="Amino Acid")

codon_table.column("rank", width=80, anchor="center")
codon_table.column("codon", width=100, anchor="center")
codon_table.column("count", width=120, anchor="center")
codon_table.column("frequency", width=120, anchor="center")
codon_table.column("amino", width=180, anchor="center")

codon_table.pack(
    side="left",
    fill="both",
    expand=True
)

scrollbar.config(
    command=codon_table.yview
)

diet_controls = tk.Frame(diet_tab)
diet_controls.pack(fill="x", padx=10, pady=10)

find_food_button = tk.Button(
    diet_controls,
    text="Find Foods",
    command=search_usda_foods,
    width=20
)

find_food_button.pack(
    side="left"
)

diet_status = tk.Label(
    diet_controls,
    text=""
)

diet_status.pack(
    side="left",
    padx=15
)

diet_frame = tk.Frame(diet_tab)

diet_frame.pack(
    fill="both",
    expand=True,
    padx=10,
    pady=5
)

diet_scrollbar = ttk.Scrollbar(
    diet_frame,
    orient="vertical"
)

diet_scrollbar.pack(
    side="right",
    fill="y"
)

diet_text = tk.Text(
    diet_frame,
    font=("Consolas", 11),
    wrap="word",
    yscrollcommand=diet_scrollbar.set
)

diet_text.pack(
    side="left",
    fill="both",
    expand=True
)

diet_scrollbar.config(
    command=diet_text.yview
)

window.mainloop()