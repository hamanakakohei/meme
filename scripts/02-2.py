#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# 注意：表示するモチーフ名を適当にいじったが、motif db次第で変えた方がいいかも

def parse_args():
    parser = argparse.ArgumentParser(description="FIMO 結果から散布図を作成")
    parser.add_argument("--fimo_result", type=Path, required=True, help="fimo.tsv")
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--ref_seq_name", type=str, required=True)
    parser.add_argument("--alt_seq_name", type=str, required=True)
    parser.add_argument("--mutation_pos", type=int, required=True)
    parser.add_argument("--overlap_only", action="store_true")
    parser.add_argument("--no_label", action="store_true")
    parser.add_argument("--highlight_threshold", type=float, default=None,
                        help="X or Y軸がこの値を超える点を赤、それ以外をグレーにする（例: 4）")
    parser.add_argument("--figsize", type=float, nargs=2, default=[12, 12])
    return parser.parse_args()


def plot_fimo_scatter(
    pivot_df: pd.DataFrame,
    out_path: Path,
    title: str = "",
    show_labels: bool = True,
    figsize=(12, 12),
    highlight_threshold: float = None
) -> None:
    plt.figure(figsize=figsize)

    # --- カラーマッピング ---
    if highlight_threshold is None:
        # デフォルト（overlap / non-overlap）
        colors = {"overlap": "red", "non-overlap": "gray"}
        for label, group in pivot_df.groupby("overlap"):
            plt.scatter(
                group["ref"], group["alt"],
                color=colors.get(label, "black"),
                label=label, alpha=0.7
            )
    else:
        # 閾値を超える点を赤、そうでない点をグレー
        high_mask = (pivot_df["ref"] >= highlight_threshold) | (pivot_df["alt"] >= highlight_threshold)
        plt.scatter(
            pivot_df.loc[~high_mask, "ref"],
            pivot_df.loc[~high_mask, "alt"],
            color="gray", alpha=0.6, label=f"< {highlight_threshold}"
        )
        plt.scatter(
            pivot_df.loc[high_mask, "ref"],
            pivot_df.loc[high_mask, "alt"],
            color="red", alpha=0.8, label=f"≥ {highlight_threshold}"
        )

    if show_labels:
        highlight = (pivot_df["ref"] >= 3.5) & ((pivot_df["ref"] - pivot_df["alt"]) >= 0.1)
        for _, row in pivot_df[highlight].iterrows():
            plt.text(row["ref"], row["alt"], row["motif_final_id"], fontsize=6, color="black")

    # --- 軸設定 ---
    plt.xlabel("-log10(p-value) 野生型")
    plt.ylabel("-log10(p-value) 変異型")

    max_val = max(pivot_df["ref"].max(), pivot_df["alt"].max())
    plt.plot([0, max_val+1], [0, max_val+1], color="black", linestyle="--", alpha=0.6)

    plt.axvline(x=4, color="gray", linestyle=":", linewidth=1, alpha=0.8)
    plt.axhline(y=4, color="gray", linestyle=":", linewidth=1, alpha=0.8)

    plt.xlim(2, max_val + 0.5)
    plt.ylim(2, max_val + 0.5)

    # 軸目盛を整数に限定
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    if title:
        plt.title(title)

    #plt.legend()
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=1000)
    plt.close()


def main():
    args = parse_args()
    out_file = args.out_dir / f"fimo_scatter.{args.figsize[0]}.{args.figsize[1]}.png"

    # 読み込み
    df = pd.read_csv(args.fimo_result, sep="\t", comment="#")

    # 前処理
    df["log10_p"] = -np.log10(df["p-value"])
    df["sequence_name"] = df["sequence_name"].replace({
        args.ref_seq_name: "ref",
        args.alt_seq_name: "alt"
    })
    df = df[df["sequence_name"].isin(["ref", "alt"])]

    # motif名をどうするか？別に他の名前にしてもよい
    df["motif_alt_id"] = df["motif_alt_id"].fillna(df["motif_id"]) # いくつかのmotif_dbではalt_idが無いので
    df["motif_final_id"] = df["motif_alt_id"].str.split(".").str[0]
    df = df.drop(columns=["motif_id", "motif_alt_id", "q-value", "matched_sequence"])

    # ピボットテーブル作成（ref/alt 並列）
    pivot_df = df.pivot_table(
        index=["motif_final_id", "start", "stop", "strand"],
        columns="sequence_name",
        values="log10_p"
    ).reset_index()

    # 欠損を 0 補完
    pivot_df = pivot_df.fillna(0)
    pivot_df.columns.name = None  # カラム名の階層をフラットに

    # オーバーラップするモチーフに限るかどうか
    pivot_df["overlap"] = pivot_df.apply(
        lambda row: "overlap" if row["start"] <= args.mutation_pos <= row["stop"] else "non-overlap",
        axis=1
    )
    if args.overlap_only:
        pivot_df = pivot_df[pivot_df["overlap"] == "overlap"]

    # 描画
    plot_fimo_scatter(
        pivot_df,
        out_file,
        show_labels=not args.no_label,
        figsize=tuple(args.figsize),
        highlight_threshold=args.highlight_threshold
    )


if __name__ == "__main__":
    main()
~
