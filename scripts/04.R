#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

# ---- 引数の定義 ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-t", "--threshold"), type = "double", default = 1e-4, help = "p-value threshold"),
  make_option(c("--seq_name"), type = "character", default = NULL),
  make_option(c("--show_labels"), action = "store_true", default = FALSE),
  make_option(c("--base_pos"), type = "integer", default = 0, help = "ゲノムポジをx軸にするための下駄"),
  make_option(c("--xlim"), type = "character", default = NULL, help = "min,max, e.g. 1000,5000"),
  make_option(c("--width"), type = "double", default = 10),
  make_option(c("--height"), type = "double", default = 4)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- main関数 ----
plot_motif_rects <- function(input_path, output_path, pval_threshold = 1e-4, seq_name = NULL, show_labels = FALSE,
                             base_pos = 0,
                             xlim = NULL, width = 10, height = 4) {
  # データ読み込み
  # p-valueフィルタ
  # -log10変換
  df <- read.table(input_path, header = TRUE, sep = "\t") %>%
    filter(p.value <= pval_threshold) %>%
    mutate(
      logp = -log10(p.value),
      start = start + base_pos,
      stop = stop + base_pos
    )
  if(!is.null(seq_name)){
    df <- df %>% filter(sequence_name == seq_name)
  }
  # --- 長方形を積み上げるためのY軸計算 ---
  df <- df %>%
    arrange(start, stop) %>%
    mutate(layer = 0)

  for (i in 2:nrow(df)) {
    overlap <- which(df$stop[1:(i - 1)] > df$start[i])
    if (length(overlap) == 0) {
      df$layer[i] <- 0
    } else {
      df$layer[i] <- max(df$layer[overlap]) + 1
    }
  }

  # --- 描画 ---
  p <- ggplot(df) +
    geom_rect(aes(
      xmin = start-1, xmax = stop,
      ymin = layer, ymax = layer + 0.8,
      fill = logp
    ), color = "black", alpha = 0.8) +
    scale_fill_gradient(low = "white", high = "darkblue",
                        name = expression(-log[10](p-value))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(x = "Genomic Position", y = "Track layer") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  # ラベルを表示する場合
  if (show_labels && "motif_alt_id" %in% names(df)) {
    p <- p + geom_text_repel(
      aes(x = stop, y = layer + 0.4, label = motif_alt_id),
      hjust = 0,
      nudge_x = 0.05 * (max(df$stop) - min(df$start)),
      size = 3,
      segment.color = NA
    )
  }

  if (!is.null(xlim)) {
    lims <- as.numeric(strsplit(xlim, ",")[[1]])
    p <- p + coord_cartesian(xlim = lims)
  }

  # 保存
  ggsave(output_path, p, width = width, height = height, units = "cm", dpi = 800)
}

# ---- 実行 ----
plot_motif_rects(
  input_path = opt$input,
  output_path = opt$output,
  pval_threshold = opt$threshold,
  seq_name = opt$seq_name,
  show_labels = opt$show_labels,
  base_pos = opt$base_pos,
  xlim = opt$xlim,
  width = opt$width,
  height = opt$height
)
