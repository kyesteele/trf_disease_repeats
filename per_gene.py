import subprocess
import tempfile
import os
import re
import matplotlib.pyplot as plt

# configs
BINARY = "./target/release/trf_disease_repeats"
FILES_TXT = "files.txt"
MAX_PERIOD = "50"
THRESHOLDS = list(range(0, 51, 5))

NUM_CASES = 10
NUM_CONTROLS = 10

# regex patterns to parse output
repeat_re = re.compile(r"# Repeats:\s+(\d+)")
percent_re = re.compile(r"% Repeats of sequence:\s+([\d.]+)%")

with open(FILES_TXT) as f:
    files = [line.strip() for line in f if line.strip()]

assert len(files) >= NUM_CASES + NUM_CONTROLS, "Need at least 20 genes"

# run TRF on ONE gene at ONE threshold
def run_trf(file_path, threshold):
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write(file_path + "\n")
        list_path = tmp.name

    result = subprocess.run(
        [
            BINARY,
            "--input", list_path,
            "--max-period", MAX_PERIOD,
            "--threshold", str(threshold),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    os.unlink(list_path)

    repeats = 0
    percent = 0.0

    for line in result.stdout.splitlines():
        m = repeat_re.search(line)
        if m:
            repeats = int(m.group(1))
        m = percent_re.search(line)
        if m:
            percent = float(m.group(1))

    return repeats, percent

# collect per-gene curves
gene_repeat_curves = {}
gene_coverage_curves = {}

for gene in files:
    print(f"Processing {gene}")
    reps = []
    covs = []

    for t in THRESHOLDS:
        r, p = run_trf(gene, t)
        reps.append(max(r, 1))   # avoid log(0)
        covs.append(max(p, 1e-3))

    gene_repeat_curves[gene] = reps
    gene_coverage_curves[gene] = covs

# plot 1
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# expansion genes
for gene in files[:NUM_CASES]:
    axes[0].plot(
        THRESHOLDS,
        gene_repeat_curves[gene],
        alpha=0.8
    )

axes[0].set_title("Repeat-expansion genes")
axes[0].set_xlabel("Threshold")
axes[0].set_ylabel("# Repeats (log scale)")
axes[0].set_yscale("log")
axes[0].grid(True)

# control genes
for gene in files[NUM_CASES:NUM_CASES + NUM_CONTROLS]:
    axes[1].plot(
        THRESHOLDS,
        gene_repeat_curves[gene],
        alpha=0.8
    )

axes[1].set_title("Control genes")
axes[1].set_xlabel("Threshold")
axes[1].set_yscale("log")
axes[1].grid(True)

fig.suptitle("Repeat Count vs Threshold (Per Gene)")
plt.tight_layout()
plt.savefig("per_gene_repeats_vs_threshold_panels.png", dpi=300)
plt.show()

# plot 2
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# expansion genes
for gene in files[:NUM_CASES]:
    axes[0].plot(
        THRESHOLDS,
        gene_coverage_curves[gene],
        alpha=0.8
    )

axes[0].set_title("Repeat-expansion genes")
axes[0].set_xlabel("Threshold")
axes[0].set_ylabel("% Repeat Coverage (log scale)")
axes[0].set_yscale("log")
axes[0].grid(True)

# control genes
for gene in files[NUM_CASES:NUM_CASES + NUM_CONTROLS]:
    axes[1].plot(
        THRESHOLDS,
        gene_coverage_curves[gene],
        alpha=0.8
    )

axes[1].set_title("Control genes")
axes[1].set_xlabel("Threshold")
axes[1].set_yscale("log")
axes[1].grid(True)

fig.suptitle("Repeat Coverage vs Threshold (Per Gene)")
plt.tight_layout()
plt.savefig("per_gene_coverage_vs_threshold_panels.png", dpi=300)
plt.show()
