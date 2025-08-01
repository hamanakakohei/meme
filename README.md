# meme

## 🧬 概要

このスクリプトは、変異情報（VCF形式）を入力として、周囲のDNA結合モチーフにおける結合強度（-log₁₀(P-value)）の変化を、参照配列（Ref）と変異配列（Alt）で比較した散布図を作成します。


---

## ⚙️ 機能

- 変異と**オーバーラップしているモチーフ**のみを図示（オプション）
- 各ドットに**モチーフIDを表示**（オプション）
- ドットにラベルを付ける基準（しきい値）はスクリプト内にハードコーディングされています（`02-2.py`を直接編集してください）


---

## 📦 必要なもの
- [MEME Suite (FIMO)](https://meme-suite.org/)

モチーフのDBは以下からダウンロードできる：
- https://jaspar.elixir.no/downloads/
- https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.25.tgz
  

