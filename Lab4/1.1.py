import gzip
import io
import os
import random
import re
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def parse_fasta_content(content_str):
    records = []
    name = None
    parts = []

    for raw_line in content_str.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts).upper()))
            name = line[1:].strip() or "sequence"
            parts = []
        elif name is not None:
            parts.append("".join(line.split()))

    if name is not None:
        records.append((name, "".join(parts).upper()))

    valid_records = [
        (rec_name, seq) for rec_name, seq in records if seq
    ]

    if not valid_records:
        raise ValueError("No valid FASTA sequences found in data.")

    return max(valid_records, key=lambda item: len(item[1])), len(valid_records)


def read_fasta_file(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return parse_fasta_content(handle.read())
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return parse_fasta_content(handle.read())


def fetch_ncbi_fasta_by_accession(accession, log_fn=None):
    clean_acc = accession.strip().upper()
    match = re.match(r"^(GC[FA])_(\d{3})(\d{3})(\d{3})(?:\.\d+)?$", clean_acc)
    if not match:
        raise ValueError(
            "Invalid accession format. Expected format like 'GCF_001589915.1' or 'GCA_001589915'."
        )

    prefix, p1, p2, p3 = match.groups()
    dir_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{p1}/{p2}/{p3}/"

    if log_fn:
        log_fn(f"Accessing NCBI FTP directory: {dir_url}")
    time.sleep(0.4)

    req = Request(dir_url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    try:
        with urlopen(req, timeout=15) as response:
            html = response.read().decode("utf-8", errors="replace")
    except (HTTPError, URLError) as err:
        raise RuntimeError(f"Could not reach NCBI directory: {err}")

    subdirs = re.findall(rf'href=["\']({clean_acc}[^"\']*/?)["\']', html)
    if not subdirs:
        subdirs = re.findall(rf'href=["\']({prefix}_{p1}{p2}{p3}[^"\']*/?)["\']', html)

    if not subdirs:
        raise RuntimeError(f"Could not locate assembly folder for {accession} at NCBI.")

    assembly_subdir = subdirs[0].strip("/")
    assembly_url = f"{dir_url}{assembly_subdir}/"

    if log_fn:
        log_fn(f"Located assembly path: {assembly_subdir}")
    time.sleep(0.4)

    req_assembly = Request(assembly_url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urlopen(req_assembly, timeout=15) as response:
        assembly_html = response.read().decode("utf-8", errors="replace")

    files = re.findall(r'href=["\']([^"\']+_genomic\.fna\.gz)["\']', assembly_html)
    genomic_files = [
        f for f in files
        if not f.startswith("cds_") and not f.startswith("rna_") and "cds_from_genomic" not in f
    ]

    if not genomic_files:
        raise RuntimeError(f"Could not find a genomic.fna.gz file in {assembly_url}")

    target_file = genomic_files[0]
    download_url = f"{assembly_url}{target_file}"

    if log_fn:
        log_fn(f"Downloading target archive: {target_file}")
    time.sleep(0.4)

    req_file = Request(download_url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    with urlopen(req_file, timeout=30) as response:
        gz_bytes = response.read()

    if log_fn:
        log_fn(f"Decompressing {len(gz_bytes):,} bytes into memory...")
    time.sleep(0.4)

    with gzip.GzipFile(fileobj=io.BytesIO(gz_bytes)) as gz:
        decompressed_text = gz.read().decode("utf-8", errors="replace")

    return parse_fasta_content(decompressed_text), target_file


def simulate_sampling(sequence, sample_count=20000, min_len=100, max_len=150, rng=None, log_fn=None):
    if rng is None:
        rng = random.Random(42)

    seq_len = len(sequence)
    samples = []
    sampled_mask = bytearray(seq_len)

    if log_fn:
        log_fn(f"Generating {sample_count:,} random sample windows ({min_len}–{max_len} nt)...")

    step = sample_count // 5
    for i in range(sample_count):
        cur_len = rng.randint(min_len, min(max_len, seq_len))
        start = rng.randint(0, seq_len - cur_len)
        sample = sequence[start:start + cur_len]
        samples.append(sample)
        sampled_mask[start:start + cur_len] = b"\x01" * cur_len

        if log_fn and (i + 1) % step == 0:
            log_fn(f"  -> Generated {i + 1:,} / {sample_count:,} samples...")
            time.sleep(0.3)

    if log_fn:
        log_fn(f"Sampling complete. Physical coverage: {sum(sampled_mask) / seq_len:.2%}")
    time.sleep(0.4)

    return samples, sum(sampled_mask)


def compute_best_overlap(seq_a, seq_b, min_overlap=15):
    max_len = min(len(seq_a), len(seq_b))
    for length in range(max_len, min_overlap - 1, -1):
        if seq_a.endswith(seq_b[:length]):
            return length
    return 0


def assemble_greedy(samples, min_overlap=15, log_fn=None):
    k = min_overlap
    kmer_to_reads = {}

    if log_fn:
        log_fn(f"Indexing {k}-mers across all {len(samples):,} reads...")
    time.sleep(0.3)

    for idx, s in enumerate(samples):
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            kmer_to_reads.setdefault(km, set()).add(idx)

    used = [False] * len(samples)
    contigs = []

    sorted_indices = sorted(range(len(samples)), key=lambda i: len(samples[i]), reverse=True)

    if log_fn:
        log_fn("Starting greedy overlap assembly...")
    time.sleep(0.3)

    total_extended = 0

    for seed_idx in sorted_indices:
        if used[seed_idx]:
            continue

        used[seed_idx] = True
        contig = samples[seed_idx]

        while True:
            tail_region = contig[-150:]
            candidates = set()
            for i in range(len(tail_region) - k + 1):
                km = tail_region[i:i + k]
                candidates.update(kmer_to_reads.get(km, ()))

            best_overlap = 0
            best_idx = -1

            for c in candidates:
                if used[c]:
                    continue
                if samples[c] in contig:
                    used[c] = True
                    continue

                ov = compute_best_overlap(contig, samples[c], min_overlap)
                if ov > best_overlap:
                    best_overlap = ov
                    best_idx = c

            if best_idx != -1 and best_overlap >= min_overlap:
                contig += samples[best_idx][best_overlap:]
                used[best_idx] = True
                total_extended += 1
                if log_fn and total_extended % 100 == 0:
                    log_fn(f"  -> Merged {total_extended:,} reads | Active contig length: {len(contig):,} nt")
                    time.sleep(0.1)
            else:
                break

        while True:
            head_region = contig[:150]
            candidates = set()
            for i in range(len(head_region) - k + 1):
                km = head_region[i:i + k]
                candidates.update(kmer_to_reads.get(km, ()))

            best_overlap = 0
            best_idx = -1

            for c in candidates:
                if used[c]:
                    continue
                if samples[c] in contig:
                    used[c] = True
                    continue

                ov = compute_best_overlap(samples[c], contig, min_overlap)
                if ov > best_overlap:
                    best_overlap = ov
                    best_idx = c

            if best_idx != -1 and best_overlap >= min_overlap:
                contig = samples[best_idx][:-best_overlap] + contig
                used[best_idx] = True
                total_extended += 1
                if log_fn and total_extended % 100 == 0:
                    log_fn(f"  -> Merged {total_extended:,} reads | Active contig length: {len(contig):,} nt")
                    time.sleep(0.1)
            else:
                break

        contigs.append(contig)

    filtered = []
    for c in sorted(contigs, key=len, reverse=True):
        if not any(c in longer for longer in filtered):
            filtered.append(c)

    if log_fn:
        log_fn(f"Assembly finished: {len(filtered):,} contig(s) created | Longest: {len(filtered[0]):,} nt")
    time.sleep(0.4)

    return filtered


def run_full_pipeline(sequence, sample_count=20000, min_len=100, max_len=150, min_overlap=15, log_fn=None):
    valid_bases = {"A", "C", "G", "T"}
    clean_seq = "".join(b for b in sequence if b in valid_bases)

    if len(clean_seq) < min_len:
        raise ValueError(f"Sequence length ({len(clean_seq)} bp) is shorter than minimum sample size ({min_len} bp).")

    samples, covered_bases = simulate_sampling(
        clean_seq,
        sample_count=sample_count,
        min_len=min_len,
        max_len=max_len,
        log_fn=log_fn,
    )

    contigs = assemble_greedy(samples, min_overlap=min_overlap, log_fn=log_fn)

    if log_fn:
        log_fn("Verifying contig fidelity against reference sequence...")
    time.sleep(0.4)

    verified_mask = bytearray(len(clean_seq))
    matched_count = 0

    for contig in contigs:
        pos = clean_seq.find(contig)
        if pos != -1:
            matched_count += 1
            verified_mask[pos:pos + len(contig)] = b"\x01" * len(contig)
        else:
            doubled_ref = clean_seq + clean_seq
            if contig in doubled_ref:
                matched_count += 1

    reconstructed_bases = sum(verified_mask)
    exact = (len(contigs) == 1 and (contigs[0] == clean_seq or contigs[0] in (clean_seq + clean_seq)))

    if log_fn:
        log_fn("Pipeline execution completed successfully.")
    time.sleep(0.5)

    return contigs, covered_bases, reconstructed_bases, matched_count, exact, clean_seq


class App:
    def __init__(self, root):
        self.root = root
        self.local_path = ""

        root.title("Genome Sampling & Greedy Assembly Simulator")
        root.geometry("680x480")
        root.minsize(580, 380)

        main_frame = ttk.Frame(root, padding=12)
        main_frame.pack(fill="both", expand=True)

        source_frame = ttk.LabelFrame(main_frame, text="Select Genome Source", padding=10)
        source_frame.pack(fill="x", pady=(0, 10))

        ttk.Label(source_frame, text="NCBI Accession:").grid(row=0, column=0, sticky="w", padx=(0, 6))
        self.acc_entry = ttk.Entry(source_frame, width=28)
        self.acc_entry.grid(row=0, column=1, sticky="w", padx=(0, 8))

        self.ncbi_button = ttk.Button(
            source_frame,
            text="Download & Run NCBI",
            command=self.run_ncbi,
        )
        self.ncbi_button.grid(row=0, column=2, sticky="w")

        ttk.Separator(source_frame, orient="horizontal").grid(row=1, column=0, columnspan=3, sticky="ew", pady=8)

        self.file_label = ttk.Label(source_frame, text="No local file selected", foreground="gray")
        self.select_button = ttk.Button(
            source_frame,
            text="Choose Local FASTA",
            command=self.select_file,
        )
        self.select_button.grid(row=2, column=0, sticky="w", padx=(0, 6))

        self.run_local_button = ttk.Button(
            source_frame,
            text="Run Local File",
            command=self.run_local,
            state="disabled",
        )
        self.run_local_button.grid(row=2, column=1, sticky="w", padx=(0, 8))
        self.file_label.grid(row=2, column=2, sticky="w")

        self.progress = ttk.Progressbar(main_frame, mode="indeterminate")
        self.progress.pack(fill="x", pady=(0, 10))

        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True)

        self.tab_progress = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_progress, text="Progress Log")

        self.tab_output = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_output, text="Results")

        self.log_output = tk.Text(
            self.tab_progress,
            wrap="word",
            state="disabled",
            borderwidth=1,
            relief="solid",
        )
        self.log_output.pack(fill="both", expand=True, padx=4, pady=4)

        self.output = tk.Text(
            self.tab_output,
            wrap="word",
            state="disabled",
            borderwidth=1,
            relief="solid",
        )
        self.output.pack(fill="both", expand=True, padx=4, pady=4)

    def select_file(self):
        path = filedialog.askopenfilename(
            filetypes=[
                ("FASTA Files", "*.fna *.fasta *.fa *.fna.gz *.fasta.gz"),
                ("All Files", "*.*"),
            ]
        )
        if path:
            self.local_path = path
            self.file_label.config(text=os.path.basename(path), foreground="black")
            self.run_local_button.config(state="normal")
            self.log_message(f"Selected local file: {os.path.basename(path)}")

    def show(self, text):
        self.output.config(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("1.0", text)
        self.output.config(state="disabled")

    def log_message(self, message):
        self.root.after(0, self._append_log, message)

    def _append_log(self, message):
        self.log_output.config(state="normal")
        self.log_output.insert("end", message + "\n")
        self.log_output.see("end")
        self.log_output.config(state="disabled")

    def clear_logs(self):
        self.log_output.config(state="normal")
        self.log_output.delete("1.0", "end")
        self.log_output.config(state="disabled")

    def set_ui_busy(self, is_busy):
        state = "disabled" if is_busy else "normal"
        self.ncbi_button.config(state=state)
        self.select_button.config(state=state)
        self.run_local_button.config(state="disabled" if is_busy or not self.local_path else "normal")

        if is_busy:
            self.progress.start(10)
        else:
            self.progress.stop()

    def run_ncbi(self):
        acc = self.acc_entry.get().strip()
        if not acc:
            messagebox.showwarning("Warning", "Please enter an NCBI accession (e.g. GCF_001589915.1).")
            return

        self.clear_logs()
        self.notebook.select(self.tab_progress)
        self.set_ui_busy(True)
        self.show("Running NCBI pipeline...")
        self.log_message(f"Starting pipeline for accession: {acc}")
        threading.Thread(target=self.process_ncbi, args=(acc,), daemon=True).start()

    def run_local(self):
        if not self.local_path:
            return

        self.clear_logs()
        self.notebook.select(self.tab_progress)
        self.set_ui_busy(True)
        self.show("Running local FASTA pipeline...")
        self.log_message(f"Starting pipeline for local file: {os.path.basename(self.local_path)}")
        threading.Thread(target=self.process_local, daemon=True).start()

    def process_ncbi(self, accession):
        try:
            (fasta_record, record_count), target_filename = fetch_ncbi_fasta_by_accession(
                accession,
                log_fn=self.log_message,
            )
            self.execute_pipeline(
                fasta_record,
                record_count,
                source_name=f"NCBI [{accession}] - {target_filename}",
            )
        except Exception as error:
            self.root.after(0, self.fail, str(error))

    def process_local(self):
        try:
            self.log_message("Parsing local FASTA file content...")
            time.sleep(0.3)
            fasta_record, record_count = read_fasta_file(self.local_path)
            self.execute_pipeline(
                fasta_record,
                record_count,
                source_name=f"Local: {os.path.basename(self.local_path)}",
            )
        except Exception as error:
            self.root.after(0, self.fail, str(error))

    def execute_pipeline(self, fasta_record, record_count, source_name):
        name, sequence = fasta_record
        contigs, covered, reconstructed, matched, exact, clean_seq = run_full_pipeline(
            sequence,
            sample_count=20000,
            min_len=100,
            max_len=150,
            min_overlap=15,
            log_fn=self.log_message,
        )

        longest = len(contigs[0]) if contigs else 0
        cov_pct = (covered / len(clean_seq)) * 100
        recon_pct = (reconstructed / len(clean_seq)) * 100

        result = (
            f"Source: {source_name}\n"
            f"Record Header: {name}\n"
            f"Genome Length: {len(clean_seq):,} nt\n"
            f"Total FASTA Records: {record_count:,}\n"
            f"Samples Generated: 20,000 (100–150 nt each)\n"
            f"Physical Coverage: {cov_pct:.2f}%\n\n"
            f"--- Assembly Results ---\n"
            f"Assembled Contigs: {len(contigs):,}\n"
            f"Longest Contig: {longest:,} nt\n"
            f"Matched Contigs: {matched:,} / {len(contigs):,}\n"
            f"Reconstructed Coverage: {recon_pct:.2f}%\n"
            f"Exact Full Assembly: {'Yes' if exact else 'No'}"
        )

        self.root.after(0, self.finish, result)

    def finish(self, result):
        self.set_ui_busy(False)
        self.show(result)

    def fail(self, error):
        self.set_ui_busy(False)
        self.log_message(f"Error encountered: {error}")
        messagebox.showerror("Error", error)


if __name__ == "__main__":
    window = tk.Tk()
    app = App(window)
    window.mainloop()