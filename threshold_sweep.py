import subprocess
import matplotlib.pyplot as plt
import tempfile
import os
import re

# config
BINARY = "./target/release/trf_disease_repeats"
MAX_PERIOD = "3"
THRESHOLDS = list(range(0, 51, 5))

FILES_TXT = "files.txt"
NUM_CASES = 10
NUM_CONTROLS = 10

# read file list
with open(FILES_TXT) as f:
    files = [line.strip() for line in f if line.strip()]

case_files = files[:NUM_CASES]
control_files = files[-NUM_CONTROLS:]

repeat_re = re.compile(r"# Repeats:\s+(\d+)")
percent_re = re.compile(r"% Repeats of sequence:\s+([\d.]+)%")

# helper to run trf algorithm on a sequence
def run_trf(file_subset, threshold):
    """Run trf on a subset of files and return
    (avg_repeats, avg_percent_coverage)
    """
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        for f in file_subset:
            tmp.write(f + "\n")
        tmp_path = tmp.name

    result = subprocess.run(
        [
            BINARY,
            "--input", tmp_path,
            "--max-period", MAX_PERIOD,
            "--threshold", str(threshold),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    os.unlink(tmp_path)

    repeats = []
    percents = []

    for line in result.stdout.splitlines():
        m1 = repeat_re.search(line)
        m2 = percent_re.search(line)
        if m1:
            repeats.append(int(m1.group(1)))
        if m2:
            percents.append(float(m2.group(1)))

    avg_repeats = sum(repeats) / len(repeats) if repeats else 0
    avg_percent = sum(percents) / len(percents) if percents else 0

    return avg_repeats, avg_percent

# sweep thresholds
case_repeats = []
case_percent = []
control_repeats = []
control_percent = []

for t in THRESHOLDS:
    print(f"Threshold {t}")

    cr, cp = run_trf(case_files, t)
    ctr, ctp = run_trf(control_files, t)

    case_repeats.append(cr)
    case_percent.append(cp)
    control_repeats.append(ctr)
    control_percent.append(ctp)

# plot 1
plt.figure(figsize=(7, 5))
plt.plot(THRESHOLDS, case_repeats, marker="o", label="Repeat-expansion genes")
plt.plot(THRESHOLDS, control_repeats, marker="o", label="Control genes")
plt.xlabel("Threshold")
plt.ylabel("Average # Repeats per Gene")
plt.title("Repeat Count vs Threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("repeats_vs_threshold.png", dpi=300)
plt.show()

# plot 2
plt.figure(figsize=(7, 5))
plt.plot(THRESHOLDS, case_percent, marker="o", label="Repeat-expansion genes")
plt.plot(THRESHOLDS, control_percent, marker="o", label="Control genes")
plt.xlabel("Threshold")
plt.ylabel("Average % Repeat Coverage")
plt.title("Repeat Coverage vs Threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("coverage_vs_threshold.png", dpi=300)
plt.show()
