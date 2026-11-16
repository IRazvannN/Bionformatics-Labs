import random
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
MIN_GENE_LEN = 10

TRANSPOSONS = (
    ("TE1", "CAGTAC"),
    ("TE2", "TGCAAT"),
    ("TE3", "GTCAGA"),
    ("TE4", "ACCGTT"),
    ("TE5", "CTAGGA"),
    ("TE6", "GGTTCA"),
)


def reverse_complement(sequence):
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


RESERVED_MOTIFS = {
    motif
    for _, left_ir in TRANSPOSONS
    for motif in (left_ir, reverse_complement(left_ir))
}


def random_dna(length):
    while True:
        seq = "".join(random.choice("ACGT") for _ in range(length))
        if all(motif not in seq for motif in RESERVED_MOTIFS):
            return seq


def simulate_dna():
    ir_map = {
        name: (left_ir, reverse_complement(left_ir))
        for name, left_ir in TRANSPOSONS
    }

    parts = []
    ground_truth = []
    cursor = 0

    def append_segment(seq):
        nonlocal cursor
        parts.append(seq)
        cursor += len(seq)

    # 1. Normal situation: TE1 and TE2 positioned one after the other
    append_segment(random_dna(random.randint(12, 20)))

    for name in ("TE1", "TE2"):
        left_ir, right_ir = ir_map[name]
        gene_len = random.randint(MIN_GENE_LEN, 16)
        start = cursor + 1
        append_segment(left_ir + random_dna(gene_len) + right_ir)
        end = cursor
        ground_truth.append(
            (name, start, end, end - start + 1, gene_len, "Normal (Sequential)")
        )
        append_segment(random_dna(random.randint(14, 22)))

    # 2. Overlapping transposomes: TE3 and TE4 partially overlap
    # Layout: IR3_L + part_a + IR4_L + shared_mid + IR3_R + part_b + IR4_R
    left3, right3 = ir_map["TE3"]
    left4, right4 = ir_map["TE4"]
    part_a = random_dna(random.randint(4, 7))
    shared_mid = random_dna(random.randint(5, 8))
    part_b = random_dna(random.randint(4, 7))

    start3 = cursor + 1
    append_segment(left3)
    append_segment(part_a)

    start4 = cursor + 1
    append_segment(left4)
    append_segment(shared_mid)
    append_segment(right3)
    end3 = cursor

    append_segment(part_b)
    append_segment(right4)
    end4 = cursor

    gene3_len = (end3 - len(right3)) - (start3 + len(left3) - 1)
    gene4_len = (end4 - len(right4)) - (start4 + len(left4) - 1)

    ground_truth.append(
        ("TE3", start3, end3, end3 - start3 + 1, gene3_len, "Overlapping (with TE4)")
    )
    ground_truth.append(
        ("TE4", start4, end4, end4 - start4 + 1, gene4_len, "Overlapping (with TE3)")
    )

    append_segment(random_dna(random.randint(14, 22)))

    # 3. Nested transposomes: TE6 completely inside TE5
    # Layout: IR5_L + outer_left + (IR6_L + gene6 + IR6_R) + outer_right + IR5_R
    left5, right5 = ir_map["TE5"]
    left6, right6 = ir_map["TE6"]
    outer_left = random_dna(random.randint(6, 9))
    gene6 = random_dna(random.randint(MIN_GENE_LEN, 15))
    outer_right = random_dna(random.randint(6, 9))

    start5 = cursor + 1
    append_segment(left5)
    append_segment(outer_left)

    start6 = cursor + 1
    append_segment(left6 + gene6 + right6)
    end6 = cursor

    append_segment(outer_right)
    append_segment(right5)
    end5 = cursor

    gene5_len = (end5 - len(right5)) - (start5 + len(left5) - 1)
    gene6_len = len(gene6)

    ground_truth.append(
        ("TE5", start5, end5, end5 - start5 + 1, gene5_len, "Nested (contains TE6)")
    )
    ground_truth.append(
        ("TE6", start6, end6, end6 - start6 + 1, gene6_len, "Nested (inside TE5)")
    )

    append_segment(random_dna(random.randint(12, 20)))

    full_seq = "".join(parts)
    kmer_counts = {}
    for i in range(len(full_seq) - 5):
        kmer = full_seq[i : i + 6]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    if any(kmer_counts.get(motif, 0) != 1 for motif in RESERVED_MOTIFS):
        return simulate_dna()

    return full_seq, ground_truth


def classify_transposons(items):
    classified = []

    for i, item in enumerate(items):
        name, start, end, length, gene_len, left_ir, right_ir = item
        relations = []

        for j, other in enumerate(items):
            if i == j:
                continue
            o_name, o_start, o_end, *_ = other

            if o_start < start and end < o_end:
                relations.append(f"Nested (inside {o_name})")
            elif start < o_start and o_end < end:
                relations.append(f"Nested (contains {o_name})")
            elif start < o_start < end < o_end or o_start < start < o_end < end:
                ov_start = max(start, o_start)
                ov_end = min(end, o_end)
                relations.append(
                    f"Overlapping (with {o_name} at {ov_start}-{ov_end})"
                )

        case_label = (
            ", ".join(relations) if relations else "Normal (Sequential)"
        )
        classified.append(
            (name, start, end, length, gene_len, left_ir, right_ir, case_label)
        )

    return classified


def detect_transposons(sequence, ir_len=6, min_gene_len=MIN_GENE_LEN):
    n = len(sequence)
    ir_to_name = {left_ir: name for name, left_ir in TRANSPOSONS}
    min_total_len = 2 * ir_len + min_gene_len
    raw_found = []

    for start in range(0, n - min_total_len + 1):
        left_ir = sequence[start : start + ir_len]
        if left_ir not in ir_to_name:
            continue

        expected_right_ir = reverse_complement(left_ir)

        for end in range(start + min_total_len, n + 1):
            right_ir = sequence[end - ir_len : end]
            if right_ir == expected_right_ir:
                name = ir_to_name[left_ir]
                total_len = end - start
                gene_len = total_len - 2 * ir_len
                raw_found.append(
                    (
                        name,
                        start + 1,
                        end,
                        total_len,
                        gene_len,
                        left_ir,
                        right_ir,
                    )
                )
                break

    return classify_transposons(raw_found)


class TransposonApp:
    def __init__(self, root):
        self.root = root
        self.sequence = ""
        self.ground_truth = []
        self.detected = []

        root.title("Transposome Simulator & Detector (Sequential, Overlapping, Nested)")
        root.geometry("920x640")
        root.minsize(740, 480)

        controls = ttk.Frame(root, padding=10)
        controls.pack(fill="x")

        ttk.Button(
            controls,
            text="1. Generate Random DNA",
            command=self.generate,
        ).pack(side="left")

        self.detect_button = ttk.Button(
            controls,
            text="2. Detect Transposomes",
            command=self.detect,
            state="disabled",
        )
        self.detect_button.pack(side="left", padx=(8, 0))

        paned = ttk.PanedWindow(root, orient="vertical")
        paned.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        canvas_frame = ttk.Frame(paned)
        self.canvas = tk.Canvas(
            canvas_frame,
            background="#151922",
            height=200,
            highlightthickness=0,
        )
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", self.draw_map)
        paned.add(canvas_frame, weight=2)

        self.output = ScrolledText(
            paned,
            wrap="word",
            state="disabled",
            font=("Consolas", 10),
        )
        paned.add(self.output, weight=3)

    def generate(self):
        self.sequence, self.ground_truth = simulate_dna()
        self.detected = []
        self.detect_button.config(state="normal")

        gt_lines = [
            f"{name:<5} | Position: {f'{start}-{end}':<9} | "
            f"Length: {length:<3} bp | Gene dist: {gene_len:<2} bp | Case: {case}"
            for name, start, end, length, gene_len, case in self.ground_truth
        ]

        text = (
            f"PART 1: SIMULATED RANDOM DNA SEQUENCE ({len(self.sequence)} bp)\n"
            f"Minimum distance between inverted repeats (gene): {MIN_GENE_LEN} letters\n"
            f"{'-' * 82}\n"
            f"{self.sequence}\n\n"
            f"Simulated Transposomes (Positions & Lengths):\n"
            f"{'-' * 82}\n"
            + "\n".join(gt_lines)
        )
        self.show(text)
        self.draw_map()

    def detect(self):
        self.detected = detect_transposons(self.sequence)

        gt_lines = [
            f"{name:<5} | Position: {f'{start}-{end}':<9} | "
            f"Length: {length:<3} bp | Gene dist: {gene_len:<2} bp | Case: {case}"
            for name, start, end, length, gene_len, case in self.ground_truth
        ]

        det_lines = [
            f"{name:<5} | Position: {f'{start}-{end}':<9} | "
            f"Length: {length:<3} bp | Gene dist: {gene_len:<2} bp | "
            f"IR: {l_ir}/{r_ir} | {case}"
            for name, start, end, length, gene_len, l_ir, r_ir, case in self.detected
        ]

        text = (
            f"PART 1: SIMULATED RANDOM DNA SEQUENCE ({len(self.sequence)} bp)\n"
            f"Minimum distance between inverted repeats (gene): {MIN_GENE_LEN} letters\n"
            f"{'-' * 82}\n"
            f"{self.sequence}\n\n"
            f"Simulated Ground Truth (Positions & Lengths):\n"
            + "\n".join(gt_lines)
            + f"\n\nPART 2: DETECTED TRANSPOSOMES ({len(self.detected)} found)\n"
            f"{'-' * 82}\n"
            + "\n".join(det_lines)
        )
        self.show(text)
        self.draw_map()

    def draw_map(self, event=None):
        self.canvas.delete("all")
        if not self.sequence:
            return

        width = max(self.canvas.winfo_width(), 1)
        height = max(self.canvas.winfo_height(), 1)
        left = 45
        right = width - 45
        seq_len = len(self.sequence)

        if right <= left or seq_len == 0:
            return

        def pos_to_x(pos):
            return left + ((pos - 1) / max(seq_len - 1, 1)) * (right - left)

        backbone_y = 38
        self.canvas.create_line(
            left, backbone_y, right, backbone_y, fill="#64748b", width=3
        )
        self.canvas.create_text(
            left, backbone_y - 16, text="1 bp", fill="#cbd5e1", anchor="w", font=("Arial", 9)
        )
        self.canvas.create_text(
            right,
            backbone_y - 16,
            text=f"{seq_len} bp",
            fill="#cbd5e1",
            anchor="e",
            font=("Arial", 9),
        )

        items = (
            [
                (name, start, end, length, case)
                for name, start, end, length, _, _, _, case in self.detected
            ]
            if self.detected
            else [
                (name, start, end, length, case)
                for name, start, end, length, _, case in self.ground_truth
            ]
        )

        colors = {
            "Normal": "#38bdf8",
            "Overlapping": "#f59e0b",
            "Nested": "#a855f7",
        }

        tracks = []
        for name, start, end, length, case in items:
            placed_track = 0
            for t_idx, last_end in enumerate(tracks):
                if start > last_end + 2:
                    placed_track = t_idx
                    tracks[t_idx] = end
                    break
            else:
                placed_track = len(tracks)
                tracks.append(end)

            x1 = pos_to_x(start)
            x2 = pos_to_x(end)
            y1 = 65 + placed_track * 44
            y2 = y1 + 22

            color = "#38bdf8"
            for key, c_val in colors.items():
                if case.startswith(key):
                    color = c_val
                    break

            self.canvas.create_rectangle(
                x1, y1, x2, y2, fill=color, outline="#f8fafc", width=1
            )
            self.canvas.create_text(
                (x1 + x2) / 2,
                (y1 + y2) / 2,
                text=f"{name}: {start}-{end} ({length} bp)",
                fill="#0f172a",
                font=("Arial", 8, "bold"),
            )

        legend_y = height - 18
        legend_items = [
            ("Normal (Sequential)", "#38bdf8"),
            ("Overlapping", "#f59e0b"),
            ("Nested", "#a855f7"),
        ]
        lx = left
        for label, col in legend_items:
            self.canvas.create_rectangle(
                lx, legend_y - 6, lx + 12, legend_y + 6, fill=col, outline=""
            )
            self.canvas.create_text(
                lx + 18,
                legend_y,
                text=label,
                fill="#cbd5e1",
                anchor="w",
                font=("Arial", 9),
            )
            lx += 165

    def show(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    TransposonApp(window)
    window.mainloop()
