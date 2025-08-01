# meme

## 🧬 概要

このスクリプトは、変異情報（VCF）を入力として、タンパク質のDNA結合モチーフの強度（MEME FIMOの-log₁₀(P-value)）の変化を、Ref配列とAlt配列で比較した散布図を作成します。


---

## ⚙️ 機能

- 変異と**オーバーラップしている結合モチーフ**のみを図示（--overlap_onlyオプション）
- 各ドットに**モチーフIDを表示**（--no_labelオプション）  
  （ドットにラベルを付ける基準（X軸、Y軸のしきい値）は02-2.py内にハードコーディングされており、直接編集する必要がある）


---

## 📦 必要なもの
- [MEME Suite (FIMO)](https://meme-suite.org/)
- モチーフのDBは以下からダウンロードできる：
  - https://jaspar.elixir.no/downloads/
  - https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.25.tgz
  

