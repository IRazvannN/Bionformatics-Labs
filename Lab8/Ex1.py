import math
import tkinter as tk
from tkinter import messagebox, ttk
from tkinter.scrolledtext import ScrolledText


DEFAULT_S1 = "ACCGTGAAGCCAATAC"
DEFAULT_S2 = "AGCGTGCAGCCAATAC"
DEFAULT_MATCH = "1"
DEFAULT_MISMATCH = "-1"
DEFAULT_GAP = "-2"

MAX_MATRIX_DRAW = 45
MAX_FILL_STEPS = 2500
MAX_ALT_ALIGNMENTS = 10
VALID_BASES = set("ACGTN")
EPS = 1e-9


def parse_score(raw_value, label):
    text = raw_value.strip()
    if not text:
        raise ValueError(f"Enter a numeric value for {label}.")
    try:
        value = float(text)
    except ValueError as err:
        raise ValueError(f"{label} must be a valid number.") from err
    if not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number.")
    if value.is_integer():
        return int(value)
    return value


def format_score(value):
    if isinstance(value, float) and abs(value - round(value)) < EPS:
        return str(int(round(value)))
    if isinstance(value, float):
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return str(value)


def clean_dna(raw_value, label):
    seq = "".join(raw_value.split()).upper()
    if not seq:
        raise ValueError(f"{label} cannot be empty.")
    invalid = sorted(set(seq) - VALID_BASES)
    if invalid:
        raise ValueError(
            f"{label} contains invalid characters ({', '.join(invalid)}). "
            "Use only A, C, G, T, or N."
        )
    return seq


def needleman_wunsch(s1, s2, match, mismatch, gap):
    m = len(s1)
    n = len(s2)

    score = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[[] for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        score[i][0] = i * gap
        trace[i][0] = ["U"]

    for j in range(1, n + 1):
        score[0][j] = j * gap
        trace[0][j] = ["L"]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            sub = match if s1[i - 1] == s2[j - 1] else mismatch
            diag = score[i - 1][j - 1] + sub
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = max(diag, up, left)
            dirs = []
            if abs(diag - best) <= EPS:
                dirs.append("D")
            if abs(up - best) <= EPS:
                dirs.append("U")
            if abs(left - best) <= EPS:
                dirs.append("L")

            score[i][j] = best
            trace[i][j] = dirs

    return score, trace


def primary_traceback(s1, s2, trace):
    i, j = len(s1), len(s2)
    path_cells = [(i, j)]
    moves = []

    while i > 0 or j > 0:
        dirs = trace[i][j]
        if "D" in dirs:
            direction = "D"
            next_i, next_j = i - 1, j - 1
            c1, c2 = s1[i - 1], s2[j - 1]
        elif "U" in dirs:
            direction = "U"
            next_i, next_j = i - 1, j
            c1, c2 = s1[i - 1], "-"
        else:
            direction = "L"
            next_i, next_j = i, j - 1
            c1, c2 = "-", s2[j - 1]

        moves.append(
            {
                "from": (i, j),
                "to": (next_i, next_j),
                "dir": direction,
                "c1": c1,
                "c2": c2,
            }
        )
        i, j = next_i, next_j
        path_cells.append((i, j))

    aligned1 = "".join(m["c1"] for m in reversed(moves))
    aligned2 = "".join(m["c2"] for m in reversed(moves))

    return path_cells, moves, aligned1, aligned2


def collect_optimal_alignments(s1, s2, trace, limit=MAX_ALT_ALIGNMENTS):
    m, n = len(s1), len(s2)
    results = []
    stack = [(m, n, [], [])]

    while stack and len(results) < limit:
        i, j, rev1, rev2 = stack.pop()
        if i == 0 and j == 0:
            results.append(("".join(reversed(rev1)), "".join(reversed(rev2))))
            continue

        for direction in reversed(trace[i][j]):
            if direction == "D":
                stack.append((i - 1, j - 1, rev1 + [s1[i - 1]], rev2 + [s2[j - 1]]))
            elif direction == "U":
                stack.append((i - 1, j, rev1 + [s1[i - 1]], rev2 + ["-"]))
            else:
                stack.append((i, j - 1, rev1 + ["-"], rev2 + [s2[j - 1]]))

    return results


def build_midline(aligned1, aligned2):
    symbols = []
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            symbols.append(" ")
        elif a == b:
            symbols.append("|")
        else:
            symbols.append(".")
    return "".join(symbols)


def alignment_stats(aligned1, aligned2):
    matches = 0
    mismatches = 0
    gaps = 0
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1
    total = len(aligned1)
    identity = (matches / total * 100) if total else 0.0
    return {
        "length": total,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity": identity,
    }


def format_alignment_blocks(aligned1, aligned2, block_size=42):
    midline = build_midline(aligned1, aligned2)
    blocks = []
    for offset in range(0, len(aligned1), block_size):
        chunk1 = aligned1[offset : offset + block_size]
        chunk_mid = midline[offset : offset + block_size]
        chunk2 = aligned2[offset : offset + block_size]
        start_col = offset + 1
        end_col = offset + len(chunk1)
        blocks.append(
            f"Pos {start_col:>3}-{end_col:<3}\n"
            f"S1: {chunk1}\n"
            f"    {chunk_mid}\n"
            f"S2: {chunk2}"
        )
    return "\n\n".join(blocks)


def build_steps(m, n, moves):
    steps = [("init", 0, 0)]

    if m * n <= MAX_FILL_STEPS:
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                steps.append(("fill", i, j))

    traceback_start_idx = len(steps)
    for move_idx, move in enumerate(moves):
        i, j = move["from"]
        steps.append(("trace", move_idx, (i, j)))

    steps.append(("result", m, n))
    return steps, traceback_start_idx


class NeedlemanWunschApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Needleman-Wunsch DNA Sequence Aligner")
        self.root.geometry("1340x780")
        self.root.minsize(1040, 660)

        self.has_aligned = False
        self.s1 = ""
        self.s2 = ""
        self.match = 1
        self.mismatch = -1
        self.gap = -2

        self.score = []
        self.trace = []
        self.path_cells = []
        self.moves = []
        self.aligned1 = ""
        self.aligned2 = ""
        self.alt_alignments = []
        self.steps = []
        self.traceback_start_idx = 0
        self.current_step_idx = 0

        self.is_playing = False
        self.play_job = None
        self.live_job = None
        self._updating_scale = False
        self._suppress_trace = False
        self.cell_boxes = []
        self.step_controls = []

        self._build_ui()
        self._show_ready_state()

    def _build_ui(self):
        top_frame = ttk.Frame(self.root, padding=(10, 8, 10, 6))
        top_frame.pack(fill="x")

        seq_frame = ttk.LabelFrame(top_frame, text="DNA Sequences & Scoring", padding=8)
        seq_frame.pack(fill="x")

        ttk.Label(seq_frame, text="S1:").grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.s1_var = tk.StringVar(value=DEFAULT_S1)
        self.s1_entry = ttk.Entry(seq_frame, textvariable=self.s1_var, font=("Consolas", 10))
        self.s1_entry.grid(row=0, column=1, columnspan=5, sticky="ew", padx=(0, 12))
        self.s1_len_label = ttk.Label(seq_frame, text=f"({len(DEFAULT_S1)} bp)", width=10)
        self.s1_len_label.grid(row=0, column=6, sticky="w")

        ttk.Label(seq_frame, text="S2:").grid(
            row=1, column=0, sticky="w", padx=(0, 6), pady=(6, 0)
        )
        self.s2_var = tk.StringVar(value=DEFAULT_S2)
        self.s2_entry = ttk.Entry(seq_frame, textvariable=self.s2_var, font=("Consolas", 10))
        self.s2_entry.grid(
            row=1, column=1, columnspan=5, sticky="ew", padx=(0, 12), pady=(6, 0)
        )
        self.s2_len_label = ttk.Label(
            seq_frame, text=f"({len(DEFAULT_S2)} bp)", width=10
        )
        self.s2_len_label.grid(row=1, column=6, sticky="w", pady=(6, 0))

        params_frame = ttk.Frame(seq_frame)
        params_frame.grid(row=2, column=0, columnspan=7, sticky="ew", pady=(8, 0))

        ttk.Label(params_frame, text="Match:").pack(side="left")
        self.match_var = tk.StringVar(value=DEFAULT_MATCH)
        self.match_entry = ttk.Entry(params_frame, textvariable=self.match_var, width=8)
        self.match_entry.pack(side="left", padx=(4, 14))

        ttk.Label(params_frame, text="Mismatch:").pack(side="left")
        self.mismatch_var = tk.StringVar(value=DEFAULT_MISMATCH)
        self.mismatch_entry = ttk.Entry(
            params_frame, textvariable=self.mismatch_var, width=8
        )
        self.mismatch_entry.pack(side="left", padx=(4, 14))

        ttk.Label(params_frame, text="Gap:").pack(side="left")
        self.gap_var = tk.StringVar(value=DEFAULT_GAP)
        self.gap_entry = ttk.Entry(params_frame, textvariable=self.gap_var, width=8)
        self.gap_entry.pack(side="left", padx=(4, 18))

        ttk.Button(
            params_frame,
            text="Align Sequences",
            command=lambda: self.run_alignment(show_errors=True),
        ).pack(side="left", padx=(0, 8))
        ttk.Button(
            params_frame, text="Reset Defaults", command=self.reset_defaults
        ).pack(side="left")

        self.score_summary_label = ttk.Label(
            params_frame, text="Not run yet", font=("Arial", 10, "bold")
        )
        self.score_summary_label.pack(side="right")

        seq_frame.columnconfigure(1, weight=1)

        for entry in (
            self.s1_entry,
            self.s2_entry,
            self.match_entry,
            self.mismatch_entry,
            self.gap_entry,
        ):
            entry.bind("<Return>", lambda _e: self.run_alignment(show_errors=True))

        for var in (
            self.s1_var,
            self.s2_var,
            self.match_var,
            self.mismatch_var,
            self.gap_var,
        ):
            var.trace_add("write", self._on_input_change)

        step_frame = ttk.LabelFrame(
            top_frame, text="Interactive Step-by-Step Explorer", padding=8
        )
        step_frame.pack(fill="x", pady=(6, 0))

        btn_row = ttk.Frame(step_frame)
        btn_row.pack(fill="x")

        btn_init = ttk.Button(
            btn_row, text="|< Init", width=8, command=self.jump_init
        )
        btn_init.pack(side="left", padx=(0, 4))

        btn_prev = ttk.Button(
            btn_row, text="< Prev", width=8, command=self.prev_step
        )
        btn_prev.pack(side="left", padx=(0, 4))

        self.play_button = ttk.Button(
            btn_row, text="Play", width=8, command=self.toggle_play
        )
        self.play_button.pack(side="left", padx=(0, 4))

        btn_next = ttk.Button(
            btn_row, text="Next >", width=8, command=self.next_step
        )
        btn_next.pack(side="left", padx=(0, 10))

        btn_tb = ttk.Button(
            btn_row, text="Go to Traceback", command=self.jump_traceback
        )
        btn_tb.pack(side="left", padx=(0, 4))

        btn_res = ttk.Button(
            btn_row, text="Full Result >|", command=self.jump_result
        )
        btn_res.pack(side="left", padx=(0, 12))

        ttk.Label(btn_row, text="Speed:").pack(side="left", padx=(4, 4))
        self.speed_var = tk.IntVar(value=120)
        self.speed_combo = ttk.Combobox(
            btn_row,
            width=8,
            state="readonly",
            values=("Fast", "Normal", "Slow"),
        )
        self.speed_combo.set("Normal")
        self.speed_combo.pack(side="left", padx=(0, 12))
        self.speed_combo.bind("<<ComboboxSelected>>", self._on_speed_change)

        self.step_var = tk.DoubleVar(value=0)
        self.step_scale = ttk.Scale(
            btn_row,
            from_=0,
            to=1,
            orient="horizontal",
            variable=self.step_var,
            command=self._on_slider_move,
        )
        self.step_scale.pack(side="left", fill="x", expand=True, padx=(0, 12))

        self.step_status_label = ttk.Label(
            btn_row, text="Ready", font=("Consolas", 10, "bold")
        )
        self.step_status_label.pack(side="right")

        self.step_controls = [
            btn_init,
            btn_prev,
            self.play_button,
            btn_next,
            btn_tb,
            btn_res,
            self.step_scale,
        ]
        self._set_step_controls_state(enabled=False)

        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        matrix_container = ttk.Frame(paned)
        legend_frame = tk.Frame(matrix_container, background="#151922", padx=10, pady=6)
        legend_frame.pack(fill="x")

        legend_items = [
            ("Current (i, j)", "#f59e0b"),
            ("Reading From", "#38bdf8"),
            ("Traceback Path", "#22c55e"),
            ("Computed", "#334155"),
            ("Click cell to inspect", None),
        ]
        for label_text, color in legend_items:
            if color:
                swatch = tk.Canvas(
                    legend_frame,
                    width=12,
                    height=12,
                    background="#151922",
                    highlightthickness=0,
                )
                swatch.create_rectangle(1, 1, 11, 11, fill=color, outline="#e2e8f0")
                swatch.pack(side="left", padx=(0, 4))
            tk.Label(
                legend_frame,
                text=label_text,
                background="#151922",
                foreground="#cbd5e1",
                font=("Arial", 9),
            ).pack(side="left", padx=(0, 12))

        canvas_grid = ttk.Frame(matrix_container)
        canvas_grid.pack(fill="both", expand=True)

        self.canvas = tk.Canvas(
            canvas_grid,
            background="#151922",
            highlightthickness=0,
        )
        self.v_scroll = ttk.Scrollbar(
            canvas_grid, orient="vertical", command=self.canvas.yview
        )
        self.h_scroll = ttk.Scrollbar(
            canvas_grid, orient="horizontal", command=self.canvas.xview
        )
        self.canvas.configure(
            yscrollcommand=self.v_scroll.set,
            xscrollcommand=self.h_scroll.set,
        )

        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.v_scroll.grid(row=0, column=1, sticky="ns")
        self.h_scroll.grid(row=1, column=0, sticky="ew")
        canvas_grid.rowconfigure(0, weight=1)
        canvas_grid.columnconfigure(0, weight=1)

        self.canvas.bind("<Button-1>", self._on_canvas_click)
        self.canvas.bind("<Configure>", self._on_canvas_resize)
        self.canvas.bind("<MouseWheel>", self._on_mousewheel)
        self.canvas.bind(
            "<Shift-MouseWheel>",
            lambda e: self.canvas.xview_scroll(-int(e.delta / 120), "units"),
        )

        paned.add(matrix_container, weight=3)

        right_container = ttk.LabelFrame(
            paned, text="Alignment & Step Details", padding=6
        )
        self.output = ScrolledText(
            right_container,
            width=54,
            wrap="word",
            state="disabled",
            font=("Consolas", 10),
        )
        self.output.pack(fill="both", expand=True)
        paned.add(right_container, weight=2)

        self.root.bind("<Left>", self._on_key_prev)
        self.root.bind("<Right>", self._on_key_next)

    def _set_step_controls_state(self, enabled):
        state = "normal" if enabled else "disabled"
        for widget in self.step_controls:
            widget.configure(state=state)

    def _show_ready_state(self):
        self.has_aligned = False
        self.score = []
        self.trace = []
        self.steps = []
        self.cell_boxes = []
        self._set_step_controls_state(enabled=False)
        self.score_summary_label.config(text="Not run yet", foreground="#475569")
        self.step_status_label.config(text="Ready")
        self._draw_ready_canvas()
        self.show_output(
            "READY TO ALIGN\n"
            "------------------------------------------------------------------\n"
            f"S1: {self.s1_var.get()}\n"
            f"S2: {self.s2_var.get()}\n\n"
            f"Scoring: Match = {self.match_var.get()}, "
            f"Mismatch = {self.mismatch_var.get()}, "
            f"Gap = {self.gap_var.get()}\n\n"
            "Modify the sequences or scores above if desired, then click\n"
            "'Align Sequences' to compute the alignment."
        )

    def _draw_ready_canvas(self):
        self.canvas.delete("all")
        self.canvas.configure(scrollregion=(0, 0, 0, 0))
        w = max(self.canvas.winfo_width(), 520)
        h = max(self.canvas.winfo_height(), 320)
        self.canvas.create_text(
            w / 2,
            h / 2,
            text="Click 'Align Sequences' to compute the Needleman-Wunsch matrix.",
            fill="#94a3b8",
            font=("Arial", 11),
        )

    def _on_canvas_resize(self, _event=None):
        if not self.has_aligned:
            self._draw_ready_canvas()

    def _on_key_prev(self, event):
        if isinstance(event.widget, (ttk.Entry, tk.Entry)):
            return
        self.prev_step()

    def _on_key_next(self, event):
        if isinstance(event.widget, (ttk.Entry, tk.Entry)):
            return
        self.next_step()

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(-int(event.delta / 120), "units")

    def _on_speed_change(self, _event=None):
        speeds = {"Fast": 35, "Normal": 120, "Slow": 380}
        self.speed_var.set(speeds.get(self.speed_combo.get(), 120))

    def reset_defaults(self):
        self._suppress_trace = True
        try:
            self.s1_var.set(DEFAULT_S1)
            self.s2_var.set(DEFAULT_S2)
            self.match_var.set(DEFAULT_MATCH)
            self.mismatch_var.set(DEFAULT_MISMATCH)
            self.gap_var.set(DEFAULT_GAP)
        finally:
            self._suppress_trace = False

        self.s1_len_label.config(text=f"({len(DEFAULT_S1)} bp)")
        self.s2_len_label.config(text=f"({len(DEFAULT_S2)} bp)")

        if self.has_aligned:
            self.run_alignment(show_errors=True)
        else:
            self._show_ready_state()

    def _on_input_change(self, *_args):
        s1_len = len("".join(self.s1_var.get().split()))
        s2_len = len("".join(self.s2_var.get().split()))
        self.s1_len_label.config(text=f"({s1_len} bp)")
        self.s2_len_label.config(text=f"({s2_len} bp)")

        if self._suppress_trace or not self.has_aligned:
            return

        if self.live_job is not None:
            self.root.after_cancel(self.live_job)
        self.live_job = self.root.after(40, self._run_live_alignment)

    def _run_live_alignment(self):
        self.live_job = None
        self.run_alignment(show_errors=False)

    def run_alignment(self, show_errors=True):
        if self.live_job is not None:
            self.root.after_cancel(self.live_job)
            self.live_job = None

        self.stop_play()
        try:
            s1 = clean_dna(self.s1_var.get(), "Sequence S1")
            s2 = clean_dna(self.s2_var.get(), "Sequence S2")
            match = parse_score(self.match_var.get(), "Match score")
            mismatch = parse_score(self.mismatch_var.get(), "Mismatch score")
            gap = parse_score(self.gap_var.get(), "Gap penalty")
        except ValueError as err:
            self.score_summary_label.config(
                text=f"Waiting for valid input: {err}",
                foreground="#b91c1c",
            )
            if show_errors:
                messagebox.showerror("Invalid Input", str(err), parent=self.root)
            return

        self.has_aligned = True
        self._set_step_controls_state(enabled=True)

        self.s1 = s1
        self.s2 = s2
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

        self.s1_len_label.config(text=f"({len(s1)} bp)")
        self.s2_len_label.config(text=f"({len(s2)} bp)")

        self.score, self.trace = needleman_wunsch(s1, s2, match, mismatch, gap)
        self.path_cells, self.moves, self.aligned1, self.aligned2 = primary_traceback(
            s1, s2, self.trace
        )
        self.alt_alignments = collect_optimal_alignments(s1, s2, self.trace)
        self.steps, self.traceback_start_idx = build_steps(
            len(s1), len(s2), self.moves
        )

        final_score = format_score(self.score[len(s1)][len(s2)])
        self.score_summary_label.config(
            text=f"Optimal Score: {final_score}",
            foreground="#0f172a",
        )

        self._updating_scale = True
        self.step_scale.config(from_=0, to=max(len(self.steps) - 1, 1))
        self._updating_scale = False

        self.go_to_step(len(self.steps) - 1)
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)

    def _on_slider_move(self, value):
        if self._updating_scale or not self.steps:
            return
        self.stop_play()
        self.go_to_step(int(round(float(value))))

    def _on_canvas_click(self, event):
        if not self.cell_boxes:
            return
        cx = self.canvas.canvasx(event.x)
        cy = self.canvas.canvasy(event.y)

        for i, j, x1, y1, x2, y2 in self.cell_boxes:
            if x1 <= cx <= x2 and y1 <= cy <= y2:
                self.stop_play()
                if i == 0 or j == 0:
                    self.go_to_step(0)
                    return

                target = ("fill", i, j)
                for idx, step in enumerate(self.steps):
                    if step == target:
                        self.go_to_step(idx)
                        return

                for idx in range(self.traceback_start_idx, len(self.steps)):
                    step = self.steps[idx]
                    if step[0] == "trace" and step[2] == (i, j):
                        self.go_to_step(idx)
                        return
                return

    def jump_init(self):
        self.stop_play()
        if self.steps:
            self.go_to_step(0)
            self.canvas.xview_moveto(0)
            self.canvas.yview_moveto(0)

    def jump_traceback(self):
        self.stop_play()
        if self.steps:
            self.go_to_step(self.traceback_start_idx)

    def jump_result(self):
        self.stop_play()
        if self.steps:
            self.go_to_step(len(self.steps) - 1)
            self.canvas.xview_moveto(0)
            self.canvas.yview_moveto(0)

    def prev_step(self):
        self.stop_play()
        if self.steps and self.current_step_idx > 0:
            self.go_to_step(self.current_step_idx - 1)

    def next_step(self):
        self.stop_play()
        if self.steps and self.current_step_idx < len(self.steps) - 1:
            self.go_to_step(self.current_step_idx + 1)

    def toggle_play(self):
        if not self.steps:
            return
        if self.is_playing:
            self.stop_play()
            return

        if self.current_step_idx >= len(self.steps) - 1:
            self.go_to_step(0)

        self.is_playing = True
        self.play_button.config(text="Pause")
        self._schedule_next()

    def stop_play(self):
        self.is_playing = False
        if self.play_job is not None:
            self.root.after_cancel(self.play_job)
            self.play_job = None
        self.play_button.config(text="Play")

    def _schedule_next(self):
        if not self.is_playing:
            return
        if self.current_step_idx >= len(self.steps) - 1:
            self.stop_play()
            return

        step = self.steps[self.current_step_idx]
        delay = self.speed_var.get()
        if step[0] == "trace":
            delay = max(delay * 2, 220)
        self.play_job = self.root.after(delay, self._advance_play)

    def _advance_play(self):
        self.play_job = None
        if not self.is_playing:
            return
        if self.current_step_idx < len(self.steps) - 1:
            self.go_to_step(self.current_step_idx + 1)
            self._schedule_next()
        else:
            self.stop_play()

    def go_to_step(self, idx):
        idx = max(0, min(idx, len(self.steps) - 1))
        self.current_step_idx = idx

        self._updating_scale = True
        self.step_var.set(idx)
        self._updating_scale = False

        self._update_status_label()
        self.draw_matrix()
        self.render_step_details()

    def _update_status_label(self):
        step = self.steps[self.current_step_idx]
        kind = step[0]
        total = len(self.steps)
        num = self.current_step_idx + 1

        if kind == "init":
            phase = "Initialization (Row 0 & Col 0)"
        elif kind == "fill":
            _, i, j = step
            phase = f"Matrix Fill | Pos ({i}, {j})"
        elif kind == "trace":
            _, move_idx, (i, j) = step
            phase = (
                f"Traceback (Going Back) | Pos ({i}, {j}) "
                f"[{move_idx + 1}/{len(self.moves)}]"
            )
        else:
            phase = "Full Result (Complete)"

        self.step_status_label.config(text=f"Step {num}/{total}  |  {phase}")

    def _trace_arrows(self, dirs):
        arrow_map = {"D": "\u2196", "U": "\u2191", "L": "\u2190"}
        return "".join(arrow_map[d] for d in dirs)

    def draw_matrix(self):
        canvas = self.canvas
        canvas.delete("all")
        self.cell_boxes = []

        if not self.score:
            return

        m, n = len(self.s1), len(self.s2)
        if m > MAX_MATRIX_DRAW or n > MAX_MATRIX_DRAW:
            canvas.configure(scrollregion=(0, 0, 720, 220))
            canvas.create_text(
                30,
                40,
                anchor="nw",
                fill="#cbd5e1",
                font=("Arial", 11, "bold"),
                text=(
                    f"DP Matrix is {m + 1} x {n + 1} (exceeds {MAX_MATRIX_DRAW} bp canvas grid limit).\n"
                    "Full alignment, score, and step-by-step traceback positions are shown on the right."
                ),
            )
            return

        step = self.steps[self.current_step_idx]
        kind = step[0]

        cell_w = 40
        cell_h = 27
        pad_x = 44
        pad_y = 42

        active_cell = None
        pred_cells = set()
        path_highlight = set()
        next_trace_cell = None
        fill_cutoff = None

        if kind == "init":
            fill_cutoff = (0, 0)
        elif kind == "fill":
            _, ci, cj = step
            active_cell = (ci, cj)
            pred_cells = {(ci - 1, cj - 1), (ci - 1, cj), (ci, cj - 1)}
            fill_cutoff = (ci, cj)
        elif kind == "trace":
            _, move_idx, (ci, cj) = step
            active_cell = (ci, cj)
            next_trace_cell = self.moves[move_idx]["to"]
            path_highlight = set(self.path_cells[: move_idx + 1])
        else:
            path_highlight = set(self.path_cells)

        total_w = pad_x + (n + 1) * cell_w + 18
        total_h = pad_y + (m + 1) * cell_h + 18
        canvas.configure(scrollregion=(0, 0, total_w, total_h))

        for j in range(n + 1):
            x1 = pad_x + j * cell_w
            cx = x1 + cell_w / 2
            label = "-" if j == 0 else self.s2[j - 1]
            is_active_col = active_cell is not None and active_cell[1] == j
            hdr_fill = "#334155" if is_active_col else "#111622"
            col_color = "#fbbf24" if is_active_col else "#e2e8f0"
            canvas.create_rectangle(
                x1, 6, x1 + cell_w, pad_y - 4, fill=hdr_fill, outline="#1e293b"
            )
            if j > 0:
                canvas.create_text(
                    cx, 14, text=str(j), fill="#94a3b8", font=("Consolas", 8)
                )
            canvas.create_text(
                cx,
                28,
                text=label,
                fill=col_color,
                font=("Consolas", 10, "bold"),
            )

        for i in range(m + 1):
            y1 = pad_y + i * cell_h
            cy = y1 + cell_h / 2
            label = "-" if i == 0 else self.s1[i - 1]
            is_active_row = active_cell is not None and active_cell[0] == i
            hdr_fill = "#334155" if is_active_row else "#111622"
            row_color = "#fbbf24" if is_active_row else "#e2e8f0"
            canvas.create_rectangle(
                6, y1, pad_x - 4, y1 + cell_h, fill=hdr_fill, outline="#1e293b"
            )
            if i > 0:
                canvas.create_text(
                    16, cy, text=str(i), fill="#94a3b8", font=("Consolas", 8)
                )
            canvas.create_text(
                31,
                cy,
                text=label,
                fill=row_color,
                font=("Consolas", 10, "bold"),
            )

        for i in range(m + 1):
            for j in range(n + 1):
                x1 = pad_x + j * cell_w
                y1 = pad_y + i * cell_h
                x2 = x1 + cell_w
                y2 = y1 + cell_h

                computed = True
                if fill_cutoff is not None and i > 0 and j > 0:
                    ci, cj = fill_cutoff
                    if (i, j) > (ci, cj):
                        computed = False

                is_active = (i, j) == active_cell
                is_pred = (i, j) in pred_cells or (i, j) == next_trace_cell
                is_path = (i, j) in path_highlight

                if is_active:
                    fill_col = "#f59e0b"
                    outline_col = "#fef3c7"
                    outline_w = 2
                    text_col = "#0f172a"
                    arrow_col = "#451a03"
                elif is_path:
                    fill_col = "#166534"
                    outline_col = "#4ade80"
                    outline_w = 2
                    text_col = "#f0fdf4"
                    arrow_col = "#86efac"
                elif is_pred:
                    fill_col = "#0c4a6e"
                    outline_col = "#38bdf8"
                    outline_w = 2
                    text_col = "#e0f2fe"
                    arrow_col = "#7dd3fc"
                elif computed:
                    fill_col = "#1e293b"
                    outline_col = "#334155"
                    outline_w = 1
                    text_col = "#e2e8f0"
                    arrow_col = "#94a3b8"
                else:
                    fill_col = "#131720"
                    outline_col = "#1e293b"
                    outline_w = 1
                    text_col = "#475569"
                    arrow_col = "#475569"

                canvas.create_rectangle(
                    x1,
                    y1,
                    x2,
                    y2,
                    fill=fill_col,
                    outline=outline_col,
                    width=outline_w,
                )
                self.cell_boxes.append((i, j, x1, y1, x2, y2))

                if computed:
                    arrows = self._trace_arrows(self.trace[i][j])
                    if arrows:
                        canvas.create_text(
                            x1 + 3,
                            y1 + 1,
                            anchor="nw",
                            text=arrows,
                            fill=arrow_col,
                            font=("Arial", 7),
                        )
                    canvas.create_text(
                        x2 - 4,
                        y2 - 2,
                        anchor="se",
                        text=format_score(self.score[i][j]),
                        fill=text_col,
                        font=("Consolas", 9, "bold"),
                    )
                else:
                    canvas.create_text(
                        (x1 + x2) / 2,
                        (y1 + y2) / 2,
                        text="\u00b7",
                        fill=text_col,
                        font=("Consolas", 10),
                    )

        if active_cell is not None:
            self._ensure_cell_visible(
                active_cell[0],
                active_cell[1],
                pad_x,
                pad_y,
                cell_w,
                cell_h,
                total_w,
                total_h,
            )

    def _ensure_cell_visible(
        self, i, j, pad_x, pad_y, cell_w, cell_h, total_w, total_h
    ):
        view_w = self.canvas.winfo_width()
        view_h = self.canvas.winfo_height()
        if view_w <= 1 or view_h <= 1:
            return

        cx1 = pad_x + j * cell_w
        cy1 = pad_y + i * cell_h
        cx2 = cx1 + cell_w
        cy2 = cy1 + cell_h

        vx1 = self.canvas.canvasx(0)
        vy1 = self.canvas.canvasy(0)
        vx2 = vx1 + view_w
        vy2 = vy1 + view_h

        if cx1 < vx1 or cx2 > vx2:
            target_x = max(0.0, min(1.0, (cx1 - view_w / 2 + cell_w / 2) / total_w))
            self.canvas.xview_moveto(target_x)
        if cy1 < vy1 or cy2 > vy2:
            target_y = max(0.0, min(1.0, (cy1 - view_h / 2 + cell_h / 2) / total_h))
            self.canvas.yview_moveto(target_y)

    def render_step_details(self):
        step = self.steps[self.current_step_idx]
        kind = step[0]
        m, n = len(self.s1), len(self.s2)

        scoring_line = (
            f"Scoring: Match = {format_score(self.match)}, "
            f"Mismatch = {format_score(self.mismatch)}, "
            f"Gap = {format_score(self.gap)}"
        )

        if kind == "init":
            row_vals = ", ".join(
                format_score(self.score[0][j]) for j in range(min(n + 1, 10))
            )
            col_vals = ", ".join(
                format_score(self.score[i][0]) for i in range(min(m + 1, 10))
            )
            suffix_r = ", ..." if n + 1 > 10 else ""
            suffix_c = ", ..." if m + 1 > 10 else ""
            text = (
                f"PHASE 1: INITIALIZATION\n"
                f"Matrix size: {m + 1} x {n + 1}\n"
                f"{scoring_line}\n"
                f"{'-' * 52}\n"
                "Base cases before comparing bases:\n"
                "  - Origin:    score[0][0] = 0\n"
                f"  - First row: score[0][j] = j * ({format_score(self.gap)})\n"
                f"               [{row_vals}{suffix_r}]\n"
                f"               (leading gaps in S1)\n"
                f"  - First col: score[i][0] = i * ({format_score(self.gap)})\n"
                f"               [{col_vals}{suffix_c}]\n"
                f"               (leading gaps in S2)"
            )
            self.show_output(text)
            return

        if kind == "fill":
            _, i, j = step
            b1 = self.s1[i - 1]
            b2 = self.s2[j - 1]
            is_match = b1 == b2
            sub = self.match if is_match else self.mismatch
            sub_kind = "Match" if is_match else "Mismatch"

            diag_prev = self.score[i - 1][j - 1]
            up_prev = self.score[i - 1][j]
            left_prev = self.score[i][j - 1]

            diag = diag_prev + sub
            up = up_prev + self.gap
            left = left_prev + self.gap
            best = self.score[i][j]
            arrows = self._trace_arrows(self.trace[i][j])

            def mark(val):
                return "  [BEST]" if abs(val - best) <= EPS else ""

            text = (
                f"PHASE 2: MATRIX FILL\n"
                f"Current Cell: (i={i}, j={j})\n"
                f"{scoring_line}\n"
                f"{'-' * 52}\n"
                f"S1[{i}] = '{b1}'  vs  S2[{j}] = '{b2}'  ({sub_kind}: {format_score(sub)})\n\n"
                f"1. Diagonal (\u2196) from ({i - 1}, {j - 1}):\n"
                f"   {format_score(diag_prev):>4} + ({format_score(sub):>3}) = {format_score(diag):>4}{mark(diag)}\n"
                f"2. Up       (\u2191) from ({i - 1}, {j}):\n"
                f"   {format_score(up_prev):>4} + ({format_score(self.gap):>3}) = {format_score(up):>4}{mark(up)}\n"
                f"3. Left     (\u2190) from ({i}, {j - 1}):\n"
                f"   {format_score(left_prev):>4} + ({format_score(self.gap):>3}) = {format_score(left):>4}{mark(left)}\n\n"
                f"Stored Score at ({i}, {j}): {format_score(best)}\n"
                f"Traceback pointer(s): {arrows}"
            )
            self.show_output(text)
            return

        if kind == "trace":
            _, move_idx, (i, j) = step
            move = self.moves[move_idx]
            next_i, next_j = move["to"]
            direction = move["dir"]
            c1, c2 = move["c1"], move["c2"]

            if direction == "D":
                if c1 == c2:
                    action = (
                        f"\u2196 Diagonal (Match: S1[{i}]='{c1}' == S2[{j}]='{c2}', "
                        f"{format_score(self.match)})"
                    )
                else:
                    action = (
                        f"\u2196 Diagonal (Mismatch: S1[{i}]='{c1}' vs S2[{j}]='{c2}', "
                        f"{format_score(self.mismatch)})"
                    )
            elif direction == "U":
                action = (
                    f"\u2191 Up (Gap in S2: S1[{i}]='{c1}' with '-', "
                    f"{format_score(self.gap)})"
                )
            else:
                action = (
                    f"\u2190 Left (Gap in S1: '-' with S2[{j}]='{c2}', "
                    f"{format_score(self.gap)})"
                )

            built_moves = self.moves[: move_idx + 1]
            partial1 = "".join(mv["c1"] for mv in reversed(built_moves))
            partial2 = "".join(mv["c2"] for mv in reversed(built_moves))
            partial_mid = build_midline(partial1, partial2)
            cursor_line = "^" + " " * max(len(partial1) - 1, 0)

            rem1 = self.s1[:next_i] or "(none)"
            rem2 = self.s2[:next_j] or "(none)"

            text = (
                f"PHASE 3: TRACEBACK (GOING BACK)\n"
                f"Step {move_idx + 1} of {len(self.moves)}\n"
                f"{'-' * 52}\n"
                f"Current Position: ({i}, {j})  [Score: {format_score(self.score[i][j])}]\n"
                f"Next Target:      ({next_i}, {next_j})  "
                f"[Score: {format_score(self.score[next_i][next_j])}]\n\n"
                f"Step taken:\n  {action}\n\n"
                f"Remaining prefixes:\n"
                f"  S1[1..{next_i}] = {rem1}\n"
                f"  S2[1..{next_j}] = {rem2}\n\n"
                f"Alignment built so far ('^' = current position):\n"
                f"    {cursor_line}\n"
                f"S1: {partial1}\n"
                f"    {partial_mid}\n"
                f"S2: {partial2}"
            )
            self.show_output(text)
            return

        stats = alignment_stats(self.aligned1, self.aligned2)
        final_score = format_score(self.score[m][n])
        blocks = format_alignment_blocks(self.aligned1, self.aligned2)

        lines = [
            "FULL RESULT: NEEDLEMAN-WUNSCH ALIGNMENT",
            scoring_line,
            "-" * 52,
            f"Optimal Score:    {final_score}",
            f"Alignment Length: {stats['length']}",
            f"Matches:          {stats['matches']}",
            f"Mismatches:       {stats['mismatches']}",
            f"Gaps:             {stats['gaps']}",
            f"Identity:         {stats['identity']:.2f}%",
            "-" * 52,
            blocks,
        ]

        if len(self.alt_alignments) > 1:
            lines.append("")
            lines.append(
                f"Tied Optimal Alignments Found: {len(self.alt_alignments)}"
                + (
                    " (showing first 10)"
                    if len(self.alt_alignments) == MAX_ALT_ALIGNMENTS
                    else ""
                )
            )
            for alt_idx, (a1, a2) in enumerate(self.alt_alignments[1:], start=2):
                lines.append(f"\nAlternative #{alt_idx}:")
                lines.append(format_alignment_blocks(a1, a2))

        self.show_output("\n".join(lines))

    def show_output(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    NeedlemanWunschApp(window)
    window.mainloop()
