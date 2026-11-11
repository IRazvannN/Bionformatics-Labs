import random
import tkinter as tk
from tkinter import ttk

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

TRANSPOSONS = [
    ("TE1", "CAGTAC", 20),
    ("TE2", "GCTAGC", 28),
    ("TE3", "TTAACC", 18),  
    ("TE4", "CCGGAA", 22),
]


def reverse_complement(sequence):
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def simulate_dna():
    host_length = random.randint(120, 220)
    dna = random_dna(host_length)
    dr_length = 5

    te_sequences = {}
    for name, left_ir, body_length in TRANSPOSONS:
        right_ir = reverse_complement(left_ir)
        te_sequences[name] = left_ir + random_dna(body_length) + right_ir

    pos1 = random.randint(10, 35)
    pos2 = random.randint(pos1 + 45, pos1 + 85)
    pos3 = random.randint(pos2 + 12, pos2 + 28)
    pos4 = random.randint(pos2 + 85, host_length + 90)

    insertions = [
        ("TE1", pos1),
        ("TE2", pos2),
        ("TE3", pos3),
        ("TE4", pos4),
    ]

    markers = list(dna)

    for name, pos in insertions:
        target_bases = [item[0] for item in markers[pos : pos + dr_length]]
        te_seq = te_sequences[name]
        tagged_te = [f"{base}:{name}" for base in te_seq]
        markers = (
            markers[: pos + dr_length]
            + tagged_te
            + target_bases
            + markers[pos + dr_length :]
        )

    final_dna = "".join(item[0] for item in markers)

    ground_truth = []
    for name, _, _ in TRANSPOSONS:
        indices = [
            idx + 1
            for idx, item in enumerate(markers)
            if item.endswith(f":{name}")
        ]
        ground_truth.append(
            (name, indices[0], indices[-1], indices[-1] - indices[0] + 1)
        )

    return final_dna, ground_truth


def detect_transposons(sequence, ir_len=6, dr_len=5, min_len=25, max_len=120):
    detected = []
    n = len(sequence)
    ir_to_name = {ir: name for name, ir, _ in TRANSPOSONS}
    expected_lengths = {
        name: len(ir) * 2 + body_len for name, ir, body_len in TRANSPOSONS
    }

    for start in range(dr_len, n - min_len - dr_len + 1):
        left_ir = sequence[start : start + ir_len]
        expected_right_ir = reverse_complement(left_ir)
        left_dr = sequence[start - dr_len : start]

        upper_bound = min(start + max_len + 1, n - dr_len + 1)

        for end in range(start + min_len, upper_bound):
            right_ir = sequence[end - ir_len : end]
            right_dr = sequence[end : end + dr_len]

            if right_ir == expected_right_ir and left_dr == right_dr:
                name = ir_to_name.get(left_ir, "Unknown")
                span = end - start
                status = "Intact"

                if name in expected_lengths and span > expected_lengths[name]:
                    diff = span - expected_lengths[name]
                    status = f"Intersected (+{diff}b inside)"

                detected.append(
                    (
                        name,
                        start + 1,
                        end,
                        span,
                        left_dr,
                        left_ir,
                        right_ir,
                        status,
                    )
                )

    return detected


class TransposonApp:
    def __init__(self, root):
        self.sequence = ""
        self.ground_truth = []

        root.title("Transposable Elements Simulator & Detector")
        root.geometry("760x540")
        root.minsize(620, 420)

        controls = ttk.Frame(root, padding=10)
        controls.pack(fill="x")

        ttk.Button(
            controls,
            text="1. Generate Artificial DNA",
            command=self.generate,
        ).pack(side="left")

        self.detect_button = ttk.Button(
            controls,
            text="2. Detect Transposons",
            command=self.detect,
            state="disabled",
        )
        self.detect_button.pack(side="left", padx=(8, 0))

        self.output = tk.Text(
            root,
            wrap="word",
            state="disabled",
            borderwidth=0,
            font=("Consolas", 10),
        )
        self.output.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10),
        )

    def generate(self):
        self.sequence, self.ground_truth = simulate_dna()
        self.detect_button.config(state="normal")

        gt_lines = [
            f"{name}: {start}-{end} (span: {span}b)"
            for name, start, end, span in self.ground_truth
        ]

        text = (
            f"PART 1: SIMULATED DNA SEQUENCE ({len(self.sequence)} bases)\n"
            f"{'-' * 68}\n"
            f"{self.sequence}\n\n"
            f"Simulated Transposons (TE3 intersects TE2):\n"
            + "\n".join(gt_lines)
        )
        self.show(text)

    def detect(self):
        found = detect_transposons(self.sequence)

        gt_lines = [
            f"{name}: {start}-{end} (span: {span}b)"
            for name, start, end, span in self.ground_truth
        ]

        det_lines = [
            f"{name:<5}  {f'{start}-{end}':<10}  Len: {span:<4}  "
            f"DR: {dr}  IR: {l_ir}/{r_ir}  [{status}]"
            for name, start, end, span, dr, l_ir, r_ir, status in found
        ]

        text = (
            f"PART 1: SIMULATED DNA SEQUENCE ({len(self.sequence)} bases)\n"
            f"{'-' * 68}\n"
            f"{self.sequence}\n\n"
            f"Simulated Ground Truth:\n"
            + "\n".join(gt_lines)
            + f"\n\nPART 2: DETECTED TRANSPOSABLE ELEMENTS ({len(found)} found)\n"
            f"{'-' * 68}\n"
            + "\n".join(det_lines)
        )
        self.show(text)

    def show(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    TransposonApp(window)
    window.mainloop()
