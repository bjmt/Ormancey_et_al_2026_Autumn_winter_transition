"""Merge coverage output in sample-table order."""
import csv
from pathlib import Path

with open("metadata/atac_samples.tsv") as f:
    samples = list(csv.DictReader(f, delimiter="\t"))
columns = []
ids = None
for sample in samples:
    path = Path("results/atac/counts") / (sample["analysis"] + ".tsv")
    with path.open() as f:
        rows = list(csv.reader(f, delimiter="\t"))
    current = [r[0] for r in rows]
    if ids is None:
        ids = current
    elif ids != current:
        raise ValueError(f"Peak order differs in {path}")
    columns.append([r[1] for r in rows])
with open("results/atac/counts.tsv", "w") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["PeakID"] + [s["analysis"] for s in samples])
    writer.writerows(zip(ids, *columns))
