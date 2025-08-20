# meme

## 🧬 概要

このスクリプトは、変異（VCF）を入力として、タンパク質のDNA結合モチーフの強度（MEME FIMOの-log₁₀(P-value)）をRefとAlt配列で比べて、散布図にする。（To do：fasta-get-markovにヒトゲノムを与えてバックグラウンドモデルを自作して使うようにする、Q値の良い解釈を思いつけばQ値で色づけたりする（Q値はインプット配列長で変わってしまう？））


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
  

