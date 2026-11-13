import math
import tkinter as tk
from collections import Counter
from itertools import combinations
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText


ENZYMES = (
    ("EcoRI", "GAATTC", 1),
    ("BamHI", "GGATCC", 1),
    ("HindIII", "AAGCTT", 1),
    ("TaqI", "TCGA", 1),
    ("HaeIII", "GGCC", 2),
)


def read_fasta(path):
    name = None
    parts = []

    with open(path, encoding="utf-8-sig") as source:
        for raw_line in source:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    raise ValueError("Select a FASTA with one sequence record.")
                name = line[1:].strip() or Path(path).name
            elif name is None:
                raise ValueError("The FASTA must start with a > header.")
            else:
                parts.append("".join(line.split()).upper())

    sequence = "".join(parts)

    if not sequence:
        raise ValueError("No sequence found.")
    if set(sequence) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("The sequence contains invalid DNA characters.")

    return name, sequence


def cut_sites(sequence, motif, offset):
    length = len(sequence)

    if length < len(motif):
        return []

    cuts = []
    start = sequence.find(motif)

    while 0 <= start < length:
        cuts.append(start + offset)
        start = sequence.find(motif, start + 1)

    return sorted(cuts)


def fragment_lengths(length, cuts):
    if not cuts:
        return [length]

    boundaries = [0, *cuts, length]

    return [
        right - left
        for left, right in zip(boundaries, boundaries[1:])
    ]


def describe(label, cuts, fragments):
    positions = ", ".join(f"{cut:,}" for cut in cuts) or "None"
    sizes = ", ".join(f"{size:,}" for size in fragments)

    return (
        f"{label}  |  {len(cuts):,} cuts  |  {len(fragments):,} fragments\n"
        f"Cut after base: {positions}\n"
        f"Fragment lengths (bp): {sizes}"
    )


class GelApp:
    def __init__(self, root):
        self.root = root
        self.name = ""
        self.sequence = ""
        self.lanes = []
        self.enzyme_vars = {}

        root.title("Restriction Digest Gel")
        root.geometry("1200x750")
        root.minsize(750, 520)

        controls = ttk.Frame(root, padding=10)
        controls.pack(fill="x")

        ttk.Button(
            controls, text="Select FASTA", command=self.select_file
        ).pack(side="left")

        self.run_button = ttk.Button(
            controls, text="Run", command=self.run, state="disabled"
        )
        self.run_button.pack(side="left", padx=(6, 16))

        ttk.Label(controls, text="Enzymes:").pack(side="left", padx=(0, 6))


        for index, (name, motif, offset) in enumerate(ENZYMES):
            var = tk.BooleanVar(value=(index < 3))
            self.enzyme_vars[name] = var
            pattern = f"{motif[:offset]}|{motif[offset:]}"
            ttk.Checkbutton(
                controls,
                text=f"{name} ({pattern})",
                variable=var,
                command=self.on_enzyme_toggle,
            ).pack(side="left", padx=4)


        paned = ttk.PanedWindow(root, orient="vertical")
        paned.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        gel_frame = ttk.Frame(paned)
        self.canvas = tk.Canvas(
            gel_frame, background="#151922", highlightthickness=0
        )
        self.h_scroll = ttk.Scrollbar(
            gel_frame, orient="horizontal", command=self.canvas.xview
        )
        self.canvas.configure(xscrollcommand=self.h_scroll.set)

        self.canvas.pack(fill="both", expand=True)
        self.h_scroll.pack(fill="x")

        self.canvas.bind("<Configure>", self.draw_gel)
        self.canvas.bind(
            "<Shift-MouseWheel>",
            lambda e: self.canvas.xview_scroll(-int(e.delta / 120), "units"),
        )

        paned.add(gel_frame, weight=3)

        self.output = ScrolledText(paned, height=11, wrap="word", state="disabled")
        paned.add(self.output, weight=1)

    def show(self, text):
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.configure(state="disabled")

    def on_enzyme_toggle(self):
        if self.sequence:
            self.run()

    def load(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.lanes = []
        self.draw_gel()
        self.show(f"{name}\nLength: {len(sequence):,} bp")
        self.run_button.configure(state="normal")
        self.run()

    def select_file(self):
        path = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fna *.fasta *.fa"), ("All files", "*.*")]
        )

        if not path:
            return

        try:
            self.load(*read_fasta(path))
        except (OSError, ValueError) as error:
            messagebox.showerror("FASTA error", str(error), parent=self.root)

    def run(self):
        if not self.sequence:
            return

        selected = [
            (name, motif, offset)
            for name, motif, offset in ENZYMES
            if self.enzyme_vars[name].get()
        ]

        if not selected:
            messagebox.showwarning(
                "No Enzymes Selected",
                "Select at least one restriction enzyme.",
                parent=self.root,
            )
            return

        length = len(self.sequence)
        self.lanes = []
        results = [
            f"{self.name}\n{length:,} bp  |  Selected enzymes: "
            f"{', '.join(name for name, _, _ in selected)}"
        ]

        enzyme_cuts = {}
        for name, motif, offset in selected:
            enzyme_cuts[name] = cut_sites(self.sequence, motif, offset)

        for r in range(1, len(selected) + 1):
            for combo in combinations(selected, r):
                if r == 1:
                    name, motif, offset = combo[0]
                    cuts = enzyme_cuts[name]
                    pattern = f"{motif[:offset]}|{motif[offset:]}"
                    text_label = f"{name} ({pattern})"
                    lane_label = name
                else:
                    names = [item[0] for item in combo]
                    combined_set = set()
                    for item_name in names:
                        combined_set.update(enzyme_cuts[item_name])
                    cuts = sorted(combined_set)
                    text_label = " + ".join(names)
                    lane_label = " +\n".join(names)

                fragments = fragment_lengths(length, cuts)
                self.lanes.append((lane_label, fragments))
                results.append(describe(text_label, cuts, fragments))

        self.show("\n\n".join(results))
        self.draw_gel()

    def draw_gel(self, event=None):
        canvas = self.canvas
        canvas.delete("all")

        if not self.lanes:
            canvas.configure(scrollregion=(0, 0, 0, 0))
            return

        view_width = canvas.winfo_width()
        height = canvas.winfo_height()
        left = 92
        min_lane_width = 85
        total_width = max(view_width, left + len(self.lanes) * min_lane_width + 24)
        right = total_width - 20
        top = 40
        bottom = height - 95

        canvas.configure(scrollregion=(0, 0, total_width, height))

        if right <= left or bottom <= top:
            return

        sizes = [size for _, fragments in self.lanes for size in fragments]
        low = math.log10(min(sizes))
        high = math.log10(max(sizes))
        span = high - low
        lane_width = (right - left) / len(self.lanes)

        def band_y(log_size):
            fraction = (high - log_size) / span if span else 0.5
            return top + fraction * (bottom - top)

        ticks = [high - span * step / 4 for step in range(5)] if span else [high]

        for tick in ticks:
            y = band_y(tick)
            canvas.create_line(left, y, right, y, fill="#303846")
            canvas.create_text(
                left - 8, y, anchor="e", fill="#c8d2df",
                text=f"{round(10 ** tick):,} bp",
            )

        font_size = 9 if len(self.lanes) <= 10 else 8

        for index, (name, fragments) in enumerate(self.lanes):
            x = left + lane_width * (index + 0.5)
            half_band = min(28, lane_width * 0.34)

            canvas.create_rectangle(
                x - 16, 13, x + 16, 23, fill="#080b11", outline=""
            )
            canvas.create_text(
                x, bottom + 10, anchor="n", justify="center",
                text=name, fill="#f0f5fa",
                font=("Arial", font_size, "bold"),
            )

            for size, count in Counter(fragments).items():
                y = band_y(math.log10(size))
                canvas.create_line(
                    x - half_band, y, x + half_band, y,
                    fill="#d9f3ff", width=min(7, 2 + count),
                )


if __name__ == "__main__":
    window = tk.Tk()
    GelApp(window)
    window.mainloop()