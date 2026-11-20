import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText


DEFAULT_SW = "30"
DEFAULT_INPUT = ">Test_Sequence_S\nCGGAACTGATCTATCTA[A(29)]CGTAGCATCTATCGATCTATCTAGCGATCTA"

SAMPLE_PROMOTERS = (
    ">Test_Sequence_S\n"
    "CGGAACTGATCTATCTA[A(29)]CGTAGCATCTATCGATCTATCTAGCGATCTA\n"
    ">Promoter_AT_Rich_TATA\n"
    "TATAAAATATTAAATATAAATTAATTTTTAAATAAATTATAAAAGGCCTATAAATAAAAT"
    "ATTAAAATATAATTTAAACCCCATATAAATAAATTAATTTTATAAAAATAAATTTAAT\n"
    ">Promoter_CG_Island\n"
    "GGCGGCGGCGCCGCGGCGCCGCGGGCGCGGCGCGCGGCGCGCGGGGCCGCGCGCGCGCCG"
    "GCGCGCCGCGGCGGCGGCCGCGCGCGGCGGCGGCGGGCGCGCCGCGCGCGGCGGCGCG\n"
    ">Promoter_Repeat_Rich\n"
    "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATC"
    "GACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGAC"
)

VALID_BASES = set("ACGTN")
PATTERN_COLORS = (
    "#38bdf8",
    "#f59e0b",
    "#22c55e",
    "#c084fc",
    "#fb7185",
    "#2dd4bf",
    "#facc15",
    "#818cf8",
)


def expand_sequence_notation(raw_seq):
    def repl(match):
        motif = match.group(1).upper()
        count = int(match.group(2))
        if count <= 0:
            raise ValueError("Repeat count inside [motif(N)] must be positive.")
        return motif * count

    expanded = re.sub(r"\[([A-Za-z]+)\((\d+)\)\]", repl, raw_seq)
    if "[" in expanded or "]" in expanded:
        raise ValueError("Invalid bracket repeat syntax. Use e.g. [A(29)].")
    return "".join(expanded.split()).upper()


def parse_fasta_records(raw_text, validate_chars=True):
    text = raw_text.strip()
    if not text:
        raise ValueError("Enter a DNA sequence or FASTA record.")

    records = []
    current_name = None
    current_parts = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None or current_parts:
                seq = "".join(current_parts)
                if validate_chars:
                    seq = expand_sequence_notation(seq)
                else:
                    seq = "".join(seq.split()).upper()
                records.append((current_name or f"Sequence {len(records) + 1}", seq))
            current_name = line[1:].strip() or f"Sequence {len(records) + 1}"
            current_parts = []
        else:
            current_parts.append(line)

    if current_name is not None or current_parts:
        seq = "".join(current_parts)
        if validate_chars:
            seq = expand_sequence_notation(seq)
        else:
            seq = "".join(seq.split()).upper()
        records.append((current_name or f"Sequence {len(records) + 1}", seq))

    cleaned = []
    for name, seq in records:
        if not seq:
            continue
        if validate_chars:
            invalid = sorted(set(seq) - VALID_BASES)
            if invalid:
                raise ValueError(
                    f"Sequence '{name}' has invalid characters ({', '.join(invalid)}). "
                    "Use only A, C, G, T, N or [A(N)] notation."
                )
        cleaned.append((name, seq))

    if not cleaned:
        raise ValueError("No non-empty FASTA sequences were found.")

    return cleaned


def compute_cg_content(sequence):
    if not sequence:
        return 0.0
    cg = sequence.count("C") + sequence.count("G")
    return (cg / len(sequence)) * 100.0


def compute_information_content(sequence):
    n = len(sequence)
    if n < 2:
        return 0.0
    total = 0.0
    for u in range(1, n):
        overlap = n - u
        matches = sum(1 for i in range(overlap) if sequence[i] == sequence[u + i])
        total += (matches / overlap) * 100.0
    return total / (n - 1)


def compute_center_of_weight(points):
    if not points:
        return 0.0, 0.0
    cx = sum(pt[0] for pt in points) / len(points)
    cy = sum(pt[1] for pt in points) / len(points)
    return cx, cy


def analyze_sequence(name, sequence, window_size=30):
    if window_size < 2:
        raise ValueError("Sliding window length must be at least 2.")
    if len(sequence) < window_size:
        raise ValueError(
            f"Sequence '{name}' length ({len(sequence)} bp) is shorter than "
            f"sliding window length ({window_size})."
        )

    windows = []
    points = []
    for start in range(0, len(sequence) - window_size + 1):
        sub = sequence[start : start + window_size]
        cg = compute_cg_content(sub)
        ic = compute_information_content(sub)
        points.append((cg, ic))
        windows.append(
            {
                "index": start + 1,
                "start": start + 1,
                "end": start + window_size,
                "cg": cg,
                "ic": ic,
            }
        )

    center_cg, center_ic = compute_center_of_weight(points)
    return {
        "name": name,
        "sequence": sequence,
        "length": len(sequence),
        "window_size": window_size,
        "whole_cg": compute_cg_content(sequence),
        "whole_ic": compute_information_content(sequence),
        "center_cg": center_cg,
        "center_ic": center_ic,
        "windows": windows,
    }


def interpolate_color(fraction):
    t = max(0.0, min(1.0, fraction))
    r = int(56 + (192 - 56) * t)
    g = int(189 + (132 - 189) * t)
    b = int(248 + (252 - 248) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


class PromoterPickerDialog(tk.Toplevel):
    def __init__(self, parent, records, title="Choose Promoter(s)", on_browse=None):
        super().__init__(parent)
        self.title(title)
        self.geometry("680x460")
        self.minsize(520, 360)
        self.transient(parent)
        self.grab_set()

        self.records = records
        self.filtered = list(enumerate(records))
        self.selected_records = None
        self.append_mode = False
        self.browse_requested = False
        self.on_browse = on_browse

        self._build_ui()
        self._apply_filter()

    def _build_ui(self):
        top = ttk.Frame(self, padding=(10, 10, 10, 6))
        top.pack(fill="x")

        ttk.Label(top, text="Filter:").pack(side="left")
        self.search_var = tk.StringVar()
        self.search_entry = ttk.Entry(top, textvariable=self.search_var)
        self.search_entry.pack(side="left", fill="x", expand=True, padx=(6, 10))
        self.search_var.trace_add("write", lambda *_a: self._apply_filter())

        self.count_label = ttk.Label(top, text="", width=22, anchor="e")
        self.count_label.pack(side="right")

        list_frame = ttk.Frame(self, padding=(10, 0, 10, 6))
        list_frame.pack(fill="both", expand=True)

        y_scroll = ttk.Scrollbar(list_frame, orient="vertical")
        y_scroll.pack(side="right", fill="y")

        self.listbox = tk.Listbox(
            list_frame,
            selectmode="extended",
            font=("Consolas", 10),
            yscrollcommand=y_scroll.set,
        )
        self.listbox.pack(side="left", fill="both", expand=True)
        y_scroll.config(command=self.listbox.yview)
        self.listbox.bind("<Double-Button-1>", lambda _e: self._confirm_selected())
        self.listbox.bind("<Return>", lambda _e: self._confirm_selected())

        opt_row = ttk.Frame(self, padding=(10, 0, 10, 6))
        opt_row.pack(fill="x")

        self.append_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            opt_row,
            text="Append to current sequence(s) instead of replacing",
            variable=self.append_var,
        ).pack(side="left")

        btn_row = ttk.Frame(self, padding=(10, 0, 10, 10))
        btn_row.pack(fill="x")

        ttk.Button(
            btn_row, text="Load Selected", command=self._confirm_selected
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            btn_row, text="Load First 5", command=self._confirm_top5
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            btn_row, text="Load All", command=self._confirm_all
        ).pack(side="left", padx=(0, 6))

        if self.on_browse is not None:
            ttk.Button(
                btn_row, text="Browse FASTA File...", command=self._trigger_browse
            ).pack(side="left", padx=(10, 0))

        ttk.Button(btn_row, text="Cancel", command=self.destroy).pack(side="right")

        self.bind("<Escape>", lambda _e: self.destroy())
        self.search_entry.focus_set()

    def _apply_filter(self):
        query = self.search_var.get().strip().lower()
        if query:
            self.filtered = [
                (idx, rec)
                for idx, rec in enumerate(self.records)
                if query in rec[0].lower()
            ]
        else:
            self.filtered = list(enumerate(self.records))

        self.listbox.delete(0, tk.END)
        for orig_idx, (name, seq) in self.filtered:
            self.listbox.insert(tk.END, f"{orig_idx + 1}. {name} ({len(seq)} bp)")

        if self.filtered:
            self.listbox.selection_set(0)
        self.count_label.config(
            text=f"Showing {len(self.filtered)} of {len(self.records)}"
        )

    def _confirm_selected(self):
        indices = self.listbox.curselection()
        if not indices:
            messagebox.showwarning(
                "No Selection", "Select at least one promoter from the list.", parent=self
            )
            return
        self.selected_records = [self.filtered[i][1] for i in indices]
        self.append_mode = self.append_var.get()
        self.destroy()

    def _confirm_top5(self):
        if not self.filtered:
            return
        self.selected_records = [rec for _, rec in self.filtered[:5]]
        self.append_mode = self.append_var.get()
        self.destroy()

    def _confirm_all(self):
        if not self.filtered:
            return
        if len(self.filtered) > 25:
            ok = messagebox.askyesno(
                "Large Selection",
                f"Load all {len(self.filtered)} promoters? Computing thousands of promoters "
                "in pure Python may take a while.",
                parent=self,
            )
            if not ok:
                return
        self.selected_records = [rec for _, rec in self.filtered]
        self.append_mode = self.append_var.get()
        self.destroy()

    def _trigger_browse(self):
        self.browse_requested = True
        self.destroy()


class PromKappaApp:
    def __init__(self, root):
        self.root = root
        self.root.title("PromKappa DNA Pattern & Center of Weight Analyzer")
        self.root.geometry("1340x820")
        self.root.minsize(1020, 680)

        self.results = []
        self.has_run = False
        self.tooltip_window = None

        self._build_ui()
        self._show_ready_state()

    def _build_ui(self):
        top_frame = ttk.Frame(self.root, padding=(10, 8, 10, 6))
        top_frame.pack(fill="x")

        input_box = ttk.LabelFrame(
            top_frame, text="DNA Sequence / Promoter FASTA Input", padding=8
        )
        input_box.pack(fill="x")

        text_row = ttk.Frame(input_box)
        text_row.pack(fill="x")

        seq_scroll = ttk.Scrollbar(text_row, orient="vertical")
        seq_scroll.pack(side="right", fill="y")

        self.seq_text = tk.Text(
            text_row,
            height=4,
            font=("Consolas", 10),
            wrap="word",
            yscrollcommand=seq_scroll.set,
        )
        self.seq_text.pack(side="left", fill="x", expand=True)
        seq_scroll.config(command=self.seq_text.yview)
        self.seq_text.insert("1.0", DEFAULT_INPUT)

        controls_row = ttk.Frame(input_box)
        controls_row.pack(fill="x", pady=(8, 0))

        ttk.Label(controls_row, text="Sliding Window (SW):").pack(side="left")
        self.sw_var = tk.StringVar(value=DEFAULT_SW)
        self.sw_entry = ttk.Entry(controls_row, textvariable=self.sw_var, width=6)
        self.sw_entry.pack(side="left", padx=(4, 14))
        self.sw_entry.bind("<Return>", lambda _e: self.run_analysis())

        ttk.Button(
            controls_row, text="Generate Pattern", command=self.run_analysis
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            controls_row, text="Load FASTA", command=self.load_fasta_file
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            controls_row, text="Load Sample Promoters", command=self.load_samples
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            controls_row, text="Copy Expanded DNA", command=self.copy_expanded_dna
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            controls_row, text="Reset Default", command=self.reset_default
        ).pack(side="left", padx=(0, 14))

        ttk.Label(controls_row, text="Chart 1 View:").pack(side="left")
        self.view_var = tk.StringVar(value="")
        self.view_combo = ttk.Combobox(
            controls_row,
            textvariable=self.view_var,
            state="disabled",
            width=22,
        )
        self.view_combo.pack(side="left", padx=(4, 10))
        self.view_combo.bind("<<ComboboxSelected>>", lambda _e: self.redraw_charts())

        self.connect_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            controls_row,
            text="Trace path",
            variable=self.connect_var,
            command=self.redraw_charts,
        ).pack(side="left")

        self.status_label = ttk.Label(
            controls_row, text="Not run yet", font=("Arial", 9, "bold")
        )
        self.status_label.pack(side="right")

        main_paned = ttk.PanedWindow(self.root, orient="vertical")
        main_paned.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        charts_frame = ttk.Frame(main_paned)
        charts_frame.columnconfigure(0, weight=1)
        charts_frame.columnconfigure(1, weight=1)
        charts_frame.rowconfigure(0, weight=1)

        chart1_box = ttk.LabelFrame(
            charts_frame,
            text="Chart 1: DNA Pattern (C+G% vs Information Content)",
            padding=4,
        )
        chart1_box.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        self.pattern_canvas = tk.Canvas(
            chart1_box, background="#151922", highlightthickness=0
        )
        self.pattern_canvas.pack(fill="both", expand=True)
        self.pattern_canvas.bind("<Configure>", lambda _e: self.draw_pattern_chart())

        chart2_box = ttk.LabelFrame(
            charts_frame,
            text="Chart 2: Center of Weight of Each Pattern",
            padding=4,
        )
        chart2_box.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        self.center_canvas = tk.Canvas(
            chart2_box, background="#151922", highlightthickness=0
        )
        self.center_canvas.pack(fill="both", expand=True)
        self.center_canvas.bind("<Configure>", lambda _e: self.draw_center_chart())

        main_paned.add(charts_frame, weight=3)

        details_box = ttk.LabelFrame(
            main_paned,
            text="Calculated C+G%, Information Content & Center of Weight",
            padding=6,
        )
        self.output = ScrolledText(
            details_box,
            height=11,
            wrap="word",
            state="disabled",
            font=("Consolas", 10),
        )
        self.output.pack(fill="both", expand=True)
        main_paned.add(details_box, weight=2)

    def _show_ready_state(self):
        self.has_run = False
        self.results = []
        self.view_combo.configure(values=(), state="disabled")
        self.view_var.set("")
        self.status_label.config(text="Not run yet", foreground="#475569")
        self.draw_pattern_chart()
        self.draw_center_chart()
        self.show_output(
            "READY TO GENERATE DNA PATTERN\n"
            "------------------------------------------------------------------------\n"
            "Default sequence loaded (with [A(29)] expanded to 29 consecutive A bases).\n"
            "Adjust the sequence, FASTA records, or Sliding Window (SW = 30), then\n"
            "click 'Generate Pattern' to compute C+G%, IC, Center of Weight, and plots."
        )

    def reset_default(self):
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert("1.0", DEFAULT_INPUT)
        self.sw_var.set(DEFAULT_SW)
        if self.has_run:
            self.run_analysis()
        else:
            self._show_ready_state()

    def _apply_chosen_records(self, chosen, append=False):
        fasta_block = "\n".join(f">{name}\n{seq}" for name, seq in chosen)
        if append:
            existing = self.seq_text.get("1.0", tk.END).strip()
            combined = f"{existing}\n{fasta_block}" if existing else fasta_block
            self.seq_text.delete("1.0", tk.END)
            self.seq_text.insert("1.0", combined)
        else:
            self.seq_text.delete("1.0", tk.END)
            self.seq_text.insert("1.0", fasta_block)
        self.run_analysis()

    def load_samples(self):
        sample_records = parse_fasta_records(SAMPLE_PROMOTERS, validate_chars=False)
        dialog = PromoterPickerDialog(
            self.root,
            sample_records,
            title="Choose Sample Promoter(s)",
            on_browse=True,
        )
        self.root.wait_window(dialog)
        if dialog.browse_requested:
            self.load_fasta_file()
            return
        if dialog.selected_records:
            self._apply_chosen_records(
                dialog.selected_records, append=dialog.append_mode
            )

    def load_fasta_file(self):
        path = filedialog.askopenfilename(
            title="Open Promoter FASTA File",
            filetypes=[
                ("FASTA / Text", "*.fasta *.fa *.fna *.txt"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                content = handle.read()
            records = parse_fasta_records(content, validate_chars=False)
        except (OSError, ValueError) as err:
            messagebox.showerror("File Error", str(err), parent=self.root)
            return

        if len(records) == 1:
            self._apply_chosen_records(records, append=False)
            return

        file_name = path.replace("\\", "/").split("/")[-1]
        dialog = PromoterPickerDialog(
            self.root,
            records,
            title=f"Choose Promoter(s) - {file_name}",
        )
        self.root.wait_window(dialog)
        if dialog.selected_records:
            self._apply_chosen_records(
                dialog.selected_records, append=dialog.append_mode
            )

    def copy_expanded_dna(self):
        try:
            records = parse_fasta_records(self.seq_text.get("1.0", tk.END))
        except ValueError as err:
            messagebox.showerror("Invalid Sequence", str(err), parent=self.root)
            return

        if len(records) == 1:
            payload = records[0][1]
        else:
            payload = "\n".join(f">{name}\n{seq}" for name, seq in records)

        self.root.clipboard_clear()
        self.root.clipboard_append(payload)
        self.status_label.config(text="Copied expanded DNA", foreground="#0f766e")

    def run_analysis(self):
        self._hide_tooltip()
        raw_text = self.seq_text.get("1.0", tk.END)
        try:
            sw = int(self.sw_var.get().strip())
        except ValueError:
            messagebox.showerror(
                "Invalid Window", "Sliding window length must be an integer.", parent=self.root
            )
            return

        try:
            records = parse_fasta_records(raw_text)
            analyzed = [analyze_sequence(name, seq, sw) for name, seq in records]
        except ValueError as err:
            messagebox.showerror("Invalid Input", str(err), parent=self.root)
            return

        prev_count = len(self.results)
        self.results = analyzed
        self.has_run = True

        names = [res["name"] for res in self.results]
        if len(names) > 1:
            options = ["All Patterns"] + names
            self.view_combo.configure(values=options, state="readonly")
            if prev_count <= 1 or self.view_var.get() not in options:
                self.view_var.set("All Patterns")
        else:
            self.view_combo.configure(values=names, state="readonly")
            self.view_var.set(names[0])

        total_windows = sum(len(r["windows"]) for r in self.results)
        self.status_label.config(
            text=f"{len(self.results)} pattern(s) | {total_windows} windows",
            foreground="#0f172a",
        )

        self.redraw_charts()
        self.render_output()

    def redraw_charts(self):
        self._hide_tooltip()
        self.draw_pattern_chart()
        self.draw_center_chart()

    def _show_tooltip(self, canvas, event, text):
        self._hide_tooltip()
        tip = tk.Toplevel(self.root)
        tip.wm_overrideredirect(True)
        tip.wm_geometry(f"+{event.x_root + 14}+{event.y_root + 10}")
        tk.Label(
            tip,
            text=text,
            background="#0f172a",
            foreground="#f8fafc",
            relief="solid",
            borderwidth=1,
            font=("Consolas", 9),
            padx=7,
            pady=4,
            justify="left",
        ).pack()
        self.tooltip_window = tip

    def _hide_tooltip(self, _event=None):
        if self.tooltip_window is not None:
            self.tooltip_window.destroy()
            self.tooltip_window = None

    def _draw_axes(self, canvas, title, x_label, y_label):
        w = max(canvas.winfo_width(), 1)
        h = max(canvas.winfo_height(), 1)
        if w < 120 or h < 120:
            return None

        left = 58
        right = w - 26
        top = 34
        bottom = h - 46

        if right <= left or bottom <= top:
            return None

        canvas.create_text(
            (left + right) / 2,
            15,
            text=title,
            fill="#e2e8f0",
            font=("Arial", 10, "bold"),
        )

        for pct in range(0, 101, 20):
            gx = left + (pct / 100.0) * (right - left)
            gy = bottom - (pct / 100.0) * (bottom - top)

            canvas.create_line(gx, top, gx, bottom, fill="#1f2937")
            canvas.create_line(left, gy, right, gy, fill="#1f2937")

            canvas.create_text(
                gx,
                bottom + 14,
                text=f"{pct}%",
                fill="#94a3b8",
                font=("Consolas", 8),
            )
            canvas.create_text(
                left - 8,
                gy,
                text=f"{pct}",
                anchor="e",
                fill="#94a3b8",
                font=("Consolas", 8),
            )

        canvas.create_line(left, bottom, right, bottom, fill="#64748b", width=2)
        canvas.create_line(left, top, left, bottom, fill="#64748b", width=2)

        canvas.create_text(
            (left + right) / 2,
            h - 12,
            text=x_label,
            fill="#cbd5e1",
            font=("Arial", 9, "bold"),
        )
        canvas.create_text(
            18,
            (top + bottom) / 2,
            text=y_label,
            fill="#cbd5e1",
            font=("Arial", 9, "bold"),
            angle=90,
        )
        return left, top, right, bottom

    def draw_pattern_chart(self):
        canvas = self.pattern_canvas
        canvas.delete("all")

        if not self.has_run or not self.results:
            w = max(canvas.winfo_width(), 320)
            h = max(canvas.winfo_height(), 220)
            canvas.create_text(
                w / 2,
                h / 2,
                text="Click 'Generate Pattern' to plot the DNA pattern.",
                fill="#94a3b8",
                font=("Arial", 10),
            )
            return

        bounds = self._draw_axes(
            canvas,
            "DNA Pattern (C+G% vs Kappa IC)",
            "C+G Content (%)",
            "Information Content (IC)",
        )
        if bounds is None:
            return
        left, top, right, bottom = bounds

        def to_xy(cg, ic):
            px = left + (max(0.0, min(100.0, cg)) / 100.0) * (right - left)
            py = bottom - (max(0.0, min(100.0, ic)) / 100.0) * (bottom - top)
            return px, py

        selected = self.view_var.get()
        if selected == "All Patterns":
            active_list = list(enumerate(self.results))
        else:
            active_list = [
                (idx, res)
                for idx, res in enumerate(self.results)
                if res["name"] == selected
            ]
            if not active_list:
                active_list = [(0, self.results[0])]

        single_mode = len(active_list) == 1

        for idx, res in active_list:
            windows = res["windows"]
            base_color = PATTERN_COLORS[idx % len(PATTERN_COLORS)]
            pts = [to_xy(w_item["cg"], w_item["ic"]) for w_item in windows]

            if self.connect_var.get() and len(pts) > 1:
                flat = [coord for pt in pts for coord in pt]
                canvas.create_line(
                    *flat,
                    fill=base_color if not single_mode else "#475569",
                    width=1,
                )

            total_pts = max(len(windows) - 1, 1)
            for p_idx, (w_item, (px, py)) in enumerate(zip(windows, pts)):
                dot_color = (
                    interpolate_color(p_idx / total_pts) if single_mode else base_color
                )
                tag = f"pat_{idx}_{p_idx}"
                r = 3.5 if single_mode else 3
                canvas.create_oval(
                    px - r,
                    py - r,
                    px + r,
                    py + r,
                    fill=dot_color,
                    outline="#0f172a",
                    width=1,
                    tags=(tag,),
                )
                tip_msg = (
                    f"{res['name']} | SW #{w_item['index']} ({w_item['start']}-{w_item['end']})\n"
                    f"C+G% = {w_item['cg']:.2f}%  |  IC = {w_item['ic']:.2f}"
                )
                canvas.tag_bind(
                    tag,
                    "<Enter>",
                    lambda e, c=canvas, msg=tip_msg: self._show_tooltip(c, e, msg),
                )
                canvas.tag_bind(tag, "<Leave>", self._hide_tooltip)

            cx, cy = to_xy(res["center_cg"], res["center_ic"])
            c_tag = f"pat_center_{idx}"
            canvas.create_line(cx - 8, cy, cx + 8, cy, fill="#fbbf24", width=2, tags=(c_tag,))
            canvas.create_line(cx, cy - 8, cx, cy + 8, fill="#fbbf24", width=2, tags=(c_tag,))
            canvas.create_oval(
                cx - 5,
                cy - 5,
                cx + 5,
                cy + 5,
                outline="#fbbf24",
                width=2,
                tags=(c_tag,),
            )
            if single_mode:
                canvas.create_text(
                    cx + 10,
                    cy - 10,
                    anchor="sw" if cy > top + 28 else "nw",
                    text=f"Center ({res['center_cg']:.2f}, {res['center_ic']:.2f})",
                    fill="#fbbf24",
                    font=("Consolas", 8, "bold"),
                )
            center_msg = (
                f"{res['name']} - Center of Weight\n"
                f"C+G% = {res['center_cg']:.2f}%  |  IC = {res['center_ic']:.2f}"
            )
            canvas.tag_bind(
                c_tag,
                "<Enter>",
                lambda e, c=canvas, msg=center_msg: self._show_tooltip(c, e, msg),
            )
            canvas.tag_bind(c_tag, "<Leave>", self._hide_tooltip)

    def draw_center_chart(self):
        canvas = self.center_canvas
        canvas.delete("all")

        if not self.has_run or not self.results:
            w = max(canvas.winfo_width(), 320)
            h = max(canvas.winfo_height(), 220)
            canvas.create_text(
                w / 2,
                h / 2,
                text="Centers of weight will be plotted here.",
                fill="#94a3b8",
                font=("Arial", 10),
            )
            return

        bounds = self._draw_axes(
            canvas,
            "Centers of Weight Across Patterns",
            "Center C+G Content (%)",
            "Center IC",
        )
        if bounds is None:
            return
        left, top, right, bottom = bounds

        def to_xy(cg, ic):
            px = left + (max(0.0, min(100.0, cg)) / 100.0) * (right - left)
            py = bottom - (max(0.0, min(100.0, ic)) / 100.0) * (bottom - top)
            return px, py

        for idx, res in enumerate(self.results):
            color = PATTERN_COLORS[idx % len(PATTERN_COLORS)]
            cx, cy = to_xy(res["center_cg"], res["center_ic"])
            tag = f"cw_{idx}"

            canvas.create_line(
                cx, bottom, cx, cy, fill="#334155", dash=(3, 3)
            )
            canvas.create_line(
                left, cy, cx, cy, fill="#334155", dash=(3, 3)
            )
            canvas.create_oval(
                cx - 6,
                cy - 6,
                cx + 6,
                cy + 6,
                fill=color,
                outline="#f8fafc",
                width=2,
                tags=(tag,),
            )

            short_title = res["name"] if len(res["name"]) <= 22 else res["name"][:19] + "..."
            label_anchor = "w" if cx < right - 195 else "e"
            dx = 10 if label_anchor == "w" else -10
            dy = -10 if cy > top + 25 else 12
            canvas.create_text(
                cx + dx,
                cy + dy,
                anchor=label_anchor,
                text=f"{short_title} ({res['center_cg']:.2f}, {res['center_ic']:.2f})",
                fill=color,
                font=("Consolas", 8, "bold"),
                tags=(tag,),
            )

            tip_msg = (
                f"{res['name']}\n"
                f"Center of Weight: C+G% = {res['center_cg']:.2f}%, IC = {res['center_ic']:.2f}\n"
                f"Length: {res['length']} bp | Windows: {len(res['windows'])}"
            )
            canvas.tag_bind(
                tag,
                "<Enter>",
                lambda e, c=canvas, msg=tip_msg: self._show_tooltip(c, e, msg),
            )
            canvas.tag_bind(tag, "<Leave>", self._hide_tooltip)

    def render_output(self):
        lines = []
        for idx, res in enumerate(self.results, start=1):
            cgs = [w["cg"] for w in res["windows"]]
            ics = [w["ic"] for w in res["windows"]]
            lines.append(
                f"[{idx}] {res['name']}  |  Length: {res['length']} bp  |  "
                f"SW = {res['window_size']}  |  Windows = {len(res['windows'])}"
            )
            lines.append(
                f"    Whole Sequence -> C+G%: {res['whole_cg']:.2f}%  |  IC: {res['whole_ic']:.2f}"
            )
            lines.append(
                f"    Center of Weight -> X (C+G%): {res['center_cg']:.2f}%  |  "
                f"Y (IC): {res['center_ic']:.2f}"
            )
            lines.append(
                f"    Window Ranges -> C+G%: [{min(cgs):.2f}%, {max(cgs):.2f}%]  |  "
                f"IC: [{min(ics):.2f}, {max(ics):.2f}]"
            )
            lines.append("    Sliding Windows:")
            for w_item in res["windows"]:
                lines.append(
                    f"      SW {w_item['index']:>3} ({w_item['start']:>4}-{w_item['end']:<4}): "
                    f"C+G% = {w_item['cg']:6.2f}%  |  IC = {w_item['ic']:6.2f}"
                )
            lines.append("-" * 76)

        self.show_output("\n".join(lines))

    def show_output(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")


if __name__ == "__main__":
    window = tk.Tk()
    PromKappaApp(window)
    window.mainloop()
