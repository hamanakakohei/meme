#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--fa", type=Path, default=Path("results/01/out.fa"))
    parser.add_argument("--motif_db", type=Path, required=True)
    parser.add_argument("--background_model", type=str, choices=["--motif--", "--uniform--", "--nrdb--"], required=True)
    parser.add_argument("--out_dir", type=Path, default=Path("results/02"))
    parser.add_argument("--log_dir", type=Path, default=Path("logs/02"))
    return parser.parse_args()


def main():
    args = parse_args()

    motif_db_name = args.motif_db.stem
    #background_label = args.background_model.strip("-")  # e.g., "--uniform--" -> "uniform"

    fimo_out = args.out_dir / f"fimo_{motif_db_name}_{args.background_model}"
    fimo_log = args.log_dir / f"fimo_{motif_db_name}_{args.background_model}.log"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fimo_log.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fimo",
        "--no-pgc",
        "--o", str(fimo_out),
        "--bfile", args.background_model,
        "--thresh", "0.1",
        str(args.motif_db),
        str(args.fa)
    ]

    with open(fimo_log, "w") as log_fh:
        subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    main()
