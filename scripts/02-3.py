#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path


def merge_ref_mut(
    df: pd.DataFrame,
    group_cols,
    seq_col,
    value_cols
):
    # 重複チェック
    duplicates = (
        df.groupby(group_cols + [seq_col])
        .size()
        .reset_index(name="count")
        .query("count > 1")
    )
    if not duplicates.empty:
        print("❌ Error: duplicate rows detected:")
        print(duplicates)
        sys.exit(1)

    # ピボット
    df_wide = df.pivot_table(
        index = group_cols,
        columns = seq_col,
        values = value_cols,
        aggfunc = "first"
    )

    # MultiIndex列を平坦化 ("score_ref" など)
    # indexを戻す
    df_wide.columns = [f"{val}_{seq}" for val, seq in df_wide.columns]
    df_wide = df_wide.reset_index()

    return df_wide


def main():
    parser = argparse.ArgumentParser(description="ref/mut行を横並びにする")
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args()

    # 読み込み
    df = pd.read_csv(args.i, sep="\t")
    df_wide = merge_ref_mut(
        df,
        group_cols = ['motif_id','motif_alt_id','start','stop','strand'],
        seq_col = 'sequence_name',
        value_cols = ['score','p-value','q-value','matched_sequence']
    )
    # 整数に戻して見た目をキレイに
    int_cols = ["start", "stop"]
    for c in int_cols:
        df_wide[c] = df_wide[c].astype("Int64")

    # 保存
    out_path = Path(args.o)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_wide.to_csv(out_path, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
