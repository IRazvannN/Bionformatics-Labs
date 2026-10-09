import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import math

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

analysis_results = []


def read_fasta():
    raw = sequence_entry.get("1.0", tk.END).strip()

    if not raw:
        return []

    lines = raw.splitlines()

    sequences = []
    current_name = None
    current_sequence = []

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_sequence:
                sequences.append(
                    (
                        current_name or f"Sequence {len(sequences) + 1}",
                        "".join(current_sequence).upper()
                    )
                )

            current_name = line[1:].strip()
            current_sequence = []
        else:
            current_sequence.append(
                line.replace(" ", "").replace("\t", "")
            )

    if current_sequence:
        sequences.append(
            (
                current_name or f"Sequence {len(sequences) + 1}",
                "".join(current_sequence).upper()
            )
        )

    if not sequences and raw:
        dna = raw.replace("\n", "").replace(" ", "").replace("\t", "").upper()
        sequences.append(("Sequence 1", dna))

    return sequences


def find_products(dna):
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

    return aug_count, products


def analyze_all():
    global analysis_results

    sequences = read_fasta()

    if not sequences:
        messagebox.showerror(
            "Invalid input",
            "Enter at least one DNA sequence."
        )
        return

    analysis_results = []

    for name, dna in sequences:
        if not dna:
            continue

        if any(base not in "ATCG" for base in dna):
            messagebox.showerror(
                "Invalid sequence",
                f"The sequence '{name}' contains characters other than A, T, C and G."
            )
            return

        aug_count, products = find_products(dna)

        analysis_results.append({
            "name": name,
            "dna": dna,
            "aug_count": aug_count,
            "products": products
        })

    if not analysis_results:
        messagebox.showerror(
            "Invalid input",
            "No valid sequences were found."
        )
        return

    show_summary()
    show_details()
    draw_chart()


def short_name(name, length=38):
    if len(name) <= length:
        return name

    return name[:length - 3] + "..."


def wrap_label(text, width=18):
    words = text.split()
    lines = []
    current = ""

    for word in words:
        if len(current) + len(word) + 1 <= width:
            if current:
                current += " " + word
            else:
                current = word
        else:
            if current:
                lines.append(current)

            current = word

    if current:
        lines.append(current)

    if not lines:
        return text

    wrapped = "\n".join(lines[:3])

    if len(lines) > 3:
        wrapped += "..."

    return wrapped


def show_tooltip(event, text):
    hide_tooltip(event)

    tooltip = tk.Toplevel(window)
    tooltip.wm_overrideredirect(True)

    x = event.x_root + 15
    y = event.y_root + 10

    tooltip.wm_geometry(f"+{x}+{y}")

    label = tk.Label(
        tooltip,
        text=text,
        background="#ffffe0",
        relief="solid",
        borderwidth=1,
        font=("Arial", 9),
        padx=7,
        pady=5,
        wraplength=600,
        justify="left"
    )
    label.pack()

    chart_canvas.tooltip_window = tooltip


def hide_tooltip(event=None):
    tooltip = getattr(chart_canvas, "tooltip_window", None)

    if tooltip:
        tooltip.destroy()
        chart_canvas.tooltip_window = None


def show_summary():
    summary_text.delete("1.0", tk.END)

    summary_text.insert(
        tk.END,
        f"Sequences analyzed: {len(analysis_results)}\n\n"
    )

    maximum = max(
        len(result["products"])
        for result in analysis_results
    )

    winners = [
        result
        for result in analysis_results
        if len(result["products"]) == maximum
    ]

    summary_text.insert(
        tk.END,
        f"Maximum possible translation products: {maximum}\n"
    )

    summary_text.insert(
        tk.END,
        "Sequence with maximum:\n"
    )

    for winner in winners:
        summary_text.insert(
            tk.END,
            f"{winner['name']}\n"
        )

    summary_text.insert(
        tk.END,
        "\n" + "=" * 80 + "\n\n"
    )

    for number, result in enumerate(analysis_results, start=1):
        summary_text.insert(
            tk.END,
            f"Sequence {number}\n"
        )

        summary_text.insert(
            tk.END,
            f"Name: {result['name']}\n"
        )

        summary_text.insert(
            tk.END,
            f"DNA length: {len(result['dna'])}\n"
        )

        summary_text.insert(
            tk.END,
            f"Total AUG codons found: {result['aug_count']}\n"
        )

        summary_text.insert(
            tk.END,
            f"Possible translation products: {len(result['products'])}\n"
        )

        summary_text.insert(
            tk.END,
            "\n" + "-" * 80 + "\n\n"
        )


def show_details():
    result_text.delete("1.0", tk.END)

    for sequence_number, result in enumerate(analysis_results, start=1):
        result_text.insert(
            tk.END,
            f"SEQUENCE {sequence_number}\n"
        )

        result_text.insert(
            tk.END,
            f"{result['name']}\n"
        )

        result_text.insert(
            tk.END,
            f"DNA length: {len(result['dna'])}\n"
        )

        result_text.insert(
            tk.END,
            f"Total AUG codons found: {result['aug_count']}\n"
        )

        result_text.insert(
            tk.END,
            f"Possible translation products: {len(result['products'])}\n"
        )

        result_text.insert(
            tk.END,
            "\n" + "=" * 100 + "\n\n"
        )

        if not result["products"]:
            result_text.insert(
                tk.END,
                "No valid AUG-to-stop translation products were found.\n\n"
            )
            continue

        for product_number, product in enumerate(result["products"], start=1):
            result_text.insert(
                tk.END,
                f"Product {product_number}\n"
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
                "\n" + "-" * 100 + "\n\n"
            )

        result_text.insert(
            tk.END,
            "\n" + "#" * 100 + "\n\n"
        )


def draw_chart():
    hide_tooltip()

    chart_canvas.delete("all")

    if not analysis_results:
        chart_canvas.configure(
            scrollregion=(0, 0, 1000, 700)
        )
        return

    counts = [
        len(result["products"])
        for result in analysis_results
    ]

    maximum = max(counts)

    left_margin = 100
    right_margin = 120
    top_margin = 80
    bottom_margin = 220
    chart_height = 420
    bar_width = 85
    gap = 100

    plot_width = len(analysis_results) * (bar_width + gap)

    canvas_width = max(
        1000,
        left_margin + plot_width + right_margin
    )

    canvas_height = (
        top_margin
        + chart_height
        + bottom_margin
    )

    chart_canvas.configure(
        scrollregion=(
            0,
            0,
            canvas_width,
            canvas_height
        )
    )

    chart_canvas.create_text(
        canvas_width / 2,
        35,
        text="Possible Translation Products by Sequence",
        font=("Arial", 18, "bold")
    )

    x_axis_y = top_margin + chart_height
    y_axis_x = left_margin

    chart_canvas.create_line(
        y_axis_x,
        top_margin,
        y_axis_x,
        x_axis_y,
        width=2
    )

    chart_canvas.create_line(
        y_axis_x,
        x_axis_y,
        canvas_width - right_margin,
        x_axis_y,
        width=2
    )

    tick_count = 5

    if maximum == 0:
        step = 1
        top_value = 5
    else:
        step = math.ceil(maximum / tick_count)
        top_value = step * tick_count

    for i in range(tick_count + 1):
        value = i * step

        y = (
            x_axis_y
            - (value / top_value) * chart_height
        )

        chart_canvas.create_line(
            y_axis_x - 5,
            y,
            y_axis_x,
            y,
            width=2
        )

        chart_canvas.create_text(
            y_axis_x - 12,
            y,
            text=str(value),
            anchor="e",
            font=("Arial", 10)
        )

        if i > 0:
            chart_canvas.create_line(
                y_axis_x,
                y,
                canvas_width - right_margin,
                y,
                fill="#d9d9d9"
            )

    chart_canvas.create_text(
        35,
        top_margin + chart_height / 2,
        text="Number\nof\nproducts",
        font=("Arial", 10, "bold"),
        justify="center"
    )

    for index, result in enumerate(analysis_results):
        count = len(result["products"])

        x0 = (
            left_margin
            + gap / 2
            + index * (bar_width + gap)
        )

        x1 = x0 + bar_width

        if top_value > 0:
            height = (
                count / top_value
            ) * chart_height
        else:
            height = 0

        y0 = x_axis_y - height
        y1 = x_axis_y

        chart_canvas.create_rectangle(
            x0,
            y0,
            x1,
            y1,
            fill="#4A90E2",
            outline=""
        )

        chart_canvas.create_text(
            (x0 + x1) / 2,
            y0 - 15,
            text=str(count),
            font=("Arial", 11, "bold")
        )

        label = wrap_label(
            short_name(result["name"], 38),
            18
        )

        tag_name = f"label_{index}"

        chart_canvas.create_text(
            (x0 + x1) / 2,
            x_axis_y + 20,
            text=label,
            anchor="n",
            font=("Arial", 9),
            justify="center",
            width=135,
            tags=(tag_name,)
        )

        chart_canvas.tag_bind(
            tag_name,
            "<Enter>",
            lambda event, full_name=result["name"]: show_tooltip(
                event,
                full_name
            )
        )

        chart_canvas.tag_bind(
            tag_name,
            "<Leave>",
            hide_tooltip
        )


def toggle_fullscreen(event=None):
    window.attributes(
        "-fullscreen",
        not window.attributes("-fullscreen")
    )


def exit_fullscreen(event=None):
    window.attributes(
        "-fullscreen",
        False
    )


window = tk.Tk()

window.title(
    "Influenza Translation Product Analyzer"
)

window.geometry("1250x900")
window.minsize(950, 700)
window.resizable(True, True)

window.bind(
    "<F11>",
    toggle_fullscreen
)

window.bind(
    "<Escape>",
    exit_fullscreen
)

main_frame = tk.Frame(window)

main_frame.pack(
    fill="both",
    expand=True,
    padx=15,
    pady=10
)

title_label = tk.Label(
    main_frame,
    text="Influenza Translation Product Analyzer",
    font=("Arial", 19, "bold")
)

title_label.pack(
    pady=8
)

instruction_label = tk.Label(
    main_frame,
    text="Paste one or more DNA FASTA sequences:"
)

instruction_label.pack()

input_frame = tk.Frame(
    main_frame
)

input_frame.pack(
    fill="x",
    pady=8
)

input_scrollbar = ttk.Scrollbar(
    input_frame,
    orient="vertical"
)

input_scrollbar.pack(
    side="right",
    fill="y"
)

sequence_entry = tk.Text(
    input_frame,
    height=9,
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
    text="Analyze All Sequences",
    command=analyze_all,
    width=25
)

analyze_button.pack(
    pady=8
)

notebook = ttk.Notebook(
    main_frame
)

notebook.pack(
    fill="both",
    expand=True,
    pady=5
)

summary_tab = tk.Frame(
    notebook
)

chart_tab = tk.Frame(
    notebook
)

details_tab = tk.Frame(
    notebook
)

notebook.add(
    summary_tab,
    text="Summary"
)

notebook.add(
    chart_tab,
    text="Chart"
)

notebook.add(
    details_tab,
    text="Detailed Products"
)

summary_scrollbar = ttk.Scrollbar(
    summary_tab,
    orient="vertical"
)

summary_scrollbar.pack(
    side="right",
    fill="y"
)

summary_text = tk.Text(
    summary_tab,
    font=("Consolas", 10),
    wrap="word",
    yscrollcommand=summary_scrollbar.set
)

summary_text.pack(
    side="left",
    fill="both",
    expand=True
)

summary_scrollbar.config(
    command=summary_text.yview
)

chart_container = tk.Frame(
    chart_tab
)

chart_container.pack(
    fill="both",
    expand=True
)

chart_y_scrollbar = ttk.Scrollbar(
    chart_container,
    orient="vertical"
)

chart_y_scrollbar.pack(
    side="right",
    fill="y"
)

chart_x_scrollbar = ttk.Scrollbar(
    chart_container,
    orient="horizontal"
)

chart_x_scrollbar.pack(
    side="bottom",
    fill="x"
)

chart_canvas = tk.Canvas(
    chart_container,
    background="white",
    yscrollcommand=chart_y_scrollbar.set,
    xscrollcommand=chart_x_scrollbar.set
)

chart_canvas.pack(
    side="left",
    fill="both",
    expand=True
)

chart_canvas.tooltip_window = None

chart_y_scrollbar.config(
    command=chart_canvas.yview
)

chart_x_scrollbar.config(
    command=chart_canvas.xview
)

chart_canvas.bind(
    "<Configure>",
    lambda event: draw_chart()
)

details_scrollbar = ttk.Scrollbar(
    details_tab,
    orient="vertical"
)

details_scrollbar.pack(
    side="right",
    fill="y"
)

result_text = tk.Text(
    details_tab,
    font=("Consolas", 10),
    wrap="word",
    yscrollcommand=details_scrollbar.set
)

result_text.pack(
    side="left",
    fill="both",
    expand=True
)

details_scrollbar.config(
    command=result_text.yview
)

window.mainloop()