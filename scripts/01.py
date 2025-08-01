#!/usr/bin/env python3

import argparse
from pathlib import Path
import pysam
import sys
import pandas as pd

# 自作モジュール
sys.path.append( '~/github/' )
from utils.vcf import read_vcf_as_df
from utils.fasta import write_seq_as_fasta


def parse_args():
    parser = argparse.ArgumentParser(description="VCFで指定したバリアントのref & alt配列をFASTAにする")
    parser.add_argument("--vcf", default="data/variants.vcf", type=Path)
    parser.add_argument("--fasta", default="data/hg38.fasta", type=Path)
    parser.add_argument("--out", default="results/01/out.fa", type=Path)
    parser.add_argument("--margin", type=int, default=50, help="変異の左右何塩基の配列をFASTAに入れるか")
    parser.add_argument("--rowlen", type=int, default=50, help="出力するFASTAの各行の塩基数")
    return parser.parse_args()


def extract_mut_sequences(
    vcf_path: Path,
    fasta_path: Path,
    output_fasta: Path,
    margin: int = 50,
    line_length: int = 50
) -> None:
    vcf_df = read_vcf_as_df(vcf_path)

    with pysam.FastaFile(fasta_path) as genome_fa:
        for _, row in vcf_df.iterrows():
            chrom = row.get('#CHROM') or row.get('CHROM')
            pos = int(row['POS'])
            ref = row['REF']
            alt = row['ALT']

            # 周辺領域を取得
            region_start = max(1, pos - margin)
            region_end = pos + margin
            seq = genome_fa.fetch(chrom, region_start - 1, region_end)

            # 変異反映
            rel_pos = pos - region_start + 1
            alt_seq = seq[:rel_pos - 1] + alt + seq[rel_pos:]

            # ヘッダーつけて保存
            base_header = f"{chrom}:{region_start}-{region_end}_{pos}"
            write_seq_as_fasta(seq, f">{base_header}ref", output_fasta, line_length, mode='a')
            write_seq_as_fasta(alt_seq, f">{base_header}alt", output_fasta, line_length, mode='a')


if __name__ == "__main__":
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    extract_mut_sequences(
        vcf_path=args.vcf,
        fasta_path=args.fasta,
        output_fasta=args.out,
        margin=args.margin,
        line_length=args.rowlen
    )
