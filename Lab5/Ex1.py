import math
import random
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


def read_fasta(path):
    name = ""
    parts = []

    with open(path, "r", encoding="utf-8", errors="replace") as file:
        for raw_line in file:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if not name:
                    name = line[1:].strip()
                elif parts:
                    break
            elif name:
                parts.append("".join(line.split()).upper())

    sequence = "".join(parts)

    if not sequence:
        raise ValueError("No FASTA sequence found")

    if any(base not in "ACGTN" for base in sequence):
        raise ValueError("Invalid DNA characters")

    if not 1000 <= len(sequence) <= 3000:
        raise ValueError(
            "Sequence length must be between 1,000 and 3,000 bases"
        )

    return name or "DNA sequence", sequence


def random_sample_size(rng, sequence_length):
    return rng.randint(100, min(3000, sequence_length))


def create_random_samples(sequence, count=10):
    rng = random.SystemRandom()
    samples = []

    for _ in range(count):
        size = random_sample_size(rng, len(sequence))
        start = rng.randint(0, len(sequence) - size)
        samples.append(sequence[start:start + size])

    return samples


class GelApp:
    def __init__(self, root):
        self.sequence = ""
        self.samples = []

        root.title("Gel Electrophoresis")
        root.geometry("850x600")
        root.minsize(700, 500)

        controls = ttk.Frame(root, padding=10)
        controls.pack(fill="x")

        ttk.Button(
            controls,
            text="Select FASTA",
            command=self.select_file,
        ).pack(side="left")

        self.run_button = ttk.Button(
            controls,
            text="Run",
            command=self.run,
            state="disabled",
        )
        self.run_button.pack(side="left", padx=(8, 0))

        self.canvas = tk.Canvas(
            root,
            background="#111318",
            highlightthickness=0,
        )
        self.canvas.pack(fill="both", expand=True, padx=10)
        self.canvas.bind("<Configure>", self.redraw)

        self.output = tk.Text(
            root,
            height=10,
            state="disabled",
            borderwidth=0,
        )
        self.output.pack(fill="x", padx=10, pady=10)

    def select_file(self):
        path = filedialog.askopenfilename(
            filetypes=[
                ("FASTA", "*.fna *.fasta *.fa"),
                ("All files", "*.*"),
            ]
        )

        if not path:
            return

        try:
            name, self.sequence = read_fasta(path)
        except Exception as error:
            messagebox.showerror("Error", str(error))
            return

        self.samples = []
        self.canvas.delete("all")
        self.run_button.config(state="normal")

        self.show(
            f"{name}\n"
            f"Length: {len(self.sequence):,} bases"
        )

    def run(self):
        self.samples = create_random_samples(self.sequence)
        self.draw_gel()

        result = "\n".join(
            f"Sample {index}: {len(sample):,} bases"
            for index, sample in enumerate(self.samples, 1)
        )

        self.show(result)

    def draw_gel(self):
        self.canvas.delete("all")

        width = max(self.canvas.winfo_width(), 1)
        height = max(self.canvas.winfo_height(), 1)

        left = 35
        right = width - 25
        top = 35
        bottom = height - 50

        self.canvas.create_rectangle(
            left,
            top,
            right,
            bottom,
            fill="#171a21",
            outline="#646b78",
        )

        lane_width = (right - left) / 10
        lengths = [len(sample) for sample in self.samples]

        log_min = math.log(100)
        log_max = math.log(3000)

        for index, length in enumerate(lengths):
            center = left + lane_width * (index + 0.5)
            well_width = min(34, lane_width * 0.55)

            self.canvas.create_rectangle(
                center - well_width / 2,
                top + 8,
                center + well_width / 2,
                top + 18,
                fill="#050608",
                outline="#050608",
            )

            migration = (
                log_max - math.log(length)
            ) / (
                log_max - log_min
            )

            y = top + 38 + migration * (bottom - top - 55)

            self.canvas.create_line(
                center - lane_width * 0.33,
                y,
                center + lane_width * 0.33,
                y,
                fill="#e7f6ff",
                width=5,
            )

            self.canvas.create_text(
                center,
                bottom + 12,
                text=f"{index + 1}",
                fill="#ffffff",
                font=("Arial", 10, "bold"),
            )

            self.canvas.create_text(
                center,
                bottom + 28,
                text=f"{length} bp",
                fill="#b9c5d6",
                font=("Arial", 9),
            )

    def redraw(self, event=None):
        if self.samples:
            self.draw_gel()

    def show(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    GelApp(window)
    window.mainloop()