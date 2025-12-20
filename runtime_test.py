import random
import subprocess
import time
import matplotlib.pyplot as plt
import tempfile
import os

# lengths to test
SIZES = [
    1_000,
    5_000,
    10_000,
    25_000,
    50_000,
    100_000,
    250_000,
    500_000,
    1_000_000,
]

BINARY = "target/release/trf_disease_repeats"
MAX_PERIOD = "50"
THRESHOLD = "50"

def gen_seq(n):
    return ''.join(random.choice("ACGT") for _ in range(n))

def write_fna(path, seq):
    with open(path, "w") as f:
        f.write(">test\n")
        f.write(seq + "\n")

sizes = []
times = []

with tempfile.TemporaryDirectory() as tmpdir:
    for n in SIZES:
        seq_path = os.path.join(tmpdir, f"seq_{n}.fna")
        list_path = os.path.join(tmpdir, "files.txt")

        write_fna(seq_path, gen_seq(n))

        with open(list_path, "w") as f:
            f.write(seq_path + "\n")

        start = time.perf_counter()

        # run Rust TRF on sequence
        subprocess.run(
            [
                BINARY,
                "--input", list_path,
                "--max-period", MAX_PERIOD,
                "--threshold", THRESHOLD,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

        elapsed = time.perf_counter() - start

        print(f"{n} bp → {elapsed:.3f} s")

        sizes.append(n)
        times.append(elapsed)

# plot on graphs
plt.figure(figsize=(7, 5))
plt.plot(sizes, times, marker="o")
plt.xlabel("Sequence length (bp)")
plt.ylabel("Runtime (seconds)")
plt.title("TRF Runtime vs Input Size")
plt.grid(True)
plt.tight_layout()
plt.savefig("runtime_vs_size.png", dpi=300)
plt.show()
