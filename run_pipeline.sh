#!/usr/bin/env bash
set -euo pipefail


source ~/miniconda3/etc/profile.d/conda.sh
conda activate meme


# 1. バリアントを与えて、ref & altの配列をFASTAにする
mkdir -p logs/01

scripts/01.py \
  --vcf data/variants.vcf \
  --fasta data/hg38.fasta \
  --out results/01/out.fa \
  --margin 50 \
  --rowlen 50 \
  > logs/01/log 2>&1


# 2-1. モチーフDBとバックモデルを指定して、先のFASTAをMEME FIMOにかける
BACKGROUND_MODELS=(--motif-- --uniform-- --nrdb--)

for motif_db in $(cat data/motif_dbs.txt); do
  for background_model in ${BACKGROUND_MODELS[@]}; do
    echo "[INFO] Running FIMO: $motif_db with $background_model"

    scripts/02-1.py \
      --fa results/01/out.fa \
      --motif_db $motif_db \
      --background_model=$background_model \
      --out_dir results/02 \
      --log_dir logs/02
  done
done


# 2-2. ref vs altで結果のスコアを比べる散布図を描く
# この例では、1つの変異を指定しているが複数指定したいときはファイルで指定するようにする
for motif_db_path in $(cat data/motif_dbs.txt); do
  motif_db=$(basename $motif_db_path)
  motif_db="${motif_db%%.*}" 

  for background_model in ${BACKGROUND_MODELS[@]}; do
    echo "[INFO] Running FIMO: $motif_db with $background_model"

    scripts/02-2.py \
      --fimo_result results/02/fimo_${motif_db}_${background_model}/fimo.tsv \
      --out_dir results/02/fimo_${motif_db}_${background_model} \
      --ref_seq_name 'chr18:63774610-63774710_63774660ref' \
      --alt_seq_name 'chr18:63774610-63774710_63774660alt' \
      --overlap_only \
      --mutation_pos 51 # fasta内の相対的ポジ
      #--no_label \
  done
done

find results/02/ -name "fimo_scatter.png" > results/02/images.txt


# 3. 散布図を結合して1つの図にする
mkdir -p results/03

scripts/03.py \
  --images results/02/images.txt \
  --out results/03/merged.png
