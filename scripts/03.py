#!/usr/bin/env python3

from pathlib import Path
import sys
import argparse

# モジュールパスの追加
sys.path.append(str(Path.home() / "github"))
from utils.others import merge_images_and_save


def parse_args():
    parser = argparse.ArgumentParser(description="複数の画像を結合して1枚にする")
    parser.add_argument("--images", type=Path, default=Path("results/02/images.txt"))
    parser.add_argument("--out", type=Path, default=Path("results/03/merged.png"))
    return parser.parse_args()


def main():
    args = parse_args()

    with args.images.open() as f:
        image_paths = [Path(line.strip()) for line in f if line.strip()]

    merge_images_and_save(image_paths, args.out, 3)


if __name__ == "__main__":
    main()

