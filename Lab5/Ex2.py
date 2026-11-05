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


def find_tandem_repeats(sequence):
    repeats = []
    n = len(sequence)

    for unit_length in range(3, 6):
        start = 0
        while start <= n - 2 * unit_length:
            unit = sequence[start:start + unit_length]

            if "N" in unit:
                start += 1
                continue

            copies = 1
            while (
                start + (copies + 1) * unit_length <= n
                and sequence[
                    start + copies * unit_length : start + (copies + 1) * unit_length
                ] == unit
            ):
                copies += 1

            total_length = unit_length * copies

            # Only accept total repetition lengths between 6 and 10 bases
            if copies >= 2 and 6 <= total_length <= 10:
                repeated_sequence = sequence[start:start + total_length]
                repeats.append(
                    (
                        start + 1,
                        start + total_length,
                        unit,
                        copies,
                        repeated_sequence,
                    )
                )
                # Advance past this repeat to prevent phase-shifted duplicates
                start += total_length
            else:
                start += 1

    return sorted(repeats, key=lambda item: (item[0], item[1], item[2]))


class RepeatApp:
    def __init__(self, root):
        self.sequence = ""

        root.title("DNA Repeat Detector")
        root.geometry("700x520")
        root.minsize(560, 400)

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

        self.output = tk.Text(
            root,
            wrap="none",
            state="disabled",
            borderwidth=0,
        )
        self.output.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10),
        )

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

        self.run_button.config(state="normal")

        self.show(
            f"{name}\n"
            f"Length: {len(self.sequence):,} bases"
        )

    def run(self):
        repeats = find_tandem_repeats(self.sequence)

        lines = [
            f"{start}-{end}   "
            f"{unit} x {copies}   "
            f"{repeated_sequence}"
            for (
                start,
                end,
                unit,
                copies,
                repeated_sequence,
            ) in repeats
        ]

        result = f"Repeats: {len(repeats)}"

        if lines:
            result += "\n" + "\n".join(lines)

        self.show(result)

    def show(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    RepeatApp(window)
    window.mainloop()