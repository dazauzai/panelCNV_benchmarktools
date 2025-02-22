library(dplyr)
library(data.table)
library(optparse)
library(fuzzyjoin)
library(readr)
option_list <- list(
  make_option(c("-s", "--stand"), type = "character", default = NULL,
              help = "Path to the standard CNV file", metavar = "character"),
  make_option(c("-t", "--true"), type = "character", default = NULL,
              help = "Path to the true positive CNV file", metavar = "character"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "Output directory for results", metavar = "character"),
  make_option(c("-f", "--func"), type = "character", default = NULL,
              help = "Function(option : event and breakpoint)", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
func <- opt$func
# 检查参数是否提供
if (is.null(opt$stand) || is.null(opt$true) || is.null(opt$out_dir)) {
  print_help(opt_parser)
  stop("Error: Missing required arguments. Please specify --stand, --true, and --out_dir.", call. = FALSE)
}
stand <- opt$stand
standard <- fread(stand, header = T) %>%
  filter(cnv != "neutral")
true <- opt$true
true_pos <- fread(true, header = T) %>%
  filter(cnv != "neutral")
out_dir <- opt$out_dir
if (func == "event") {
  rm_header <- function(input_file, output_file = tempfile(fileext = ".bed")) {
    cmd <- sprintf("tail -n +2 %s > %s", input_file, output_file)
    system(cmd)
    return(output_file)
  }
  true_col <- colnames(true_pos)
  test_col <- colnames(standard)
  temp_rm_header_true <- rm_header(true)
  temp_rm_header_test <- rm_header(stand)
  temp_tp_fp <- tempfile(fileext = ".bed")
  temp_fp <- tempfile(fileext = ".bed")
  temp_fn <- tempfile(fileext = ".bed")
  cmd_tp_fp <- sprintf("bedtools intersect -wo -a %s -b %s > %s", temp_rm_header_test, temp_rm_header_true, temp_tp_fp)
  cmd_fn <- sprintf("bedtools intersect -v -b %s -a %s > %s", temp_rm_header_test, temp_rm_header_true, temp_fn)
  cmd_fp <- sprintf("bedtools intersect -v -a %s -b %s > %s", temp_rm_header_test, temp_rm_header_true, temp_fp)
  system(cmd_tp_fp)
  system(cmd_fp)
  system(cmd_fn)
  tp_fp <- fread(temp_tp_fp, header = F)
  tp_fp_header <- c(test_col, true_col, "intersect")
  tp_fp <- setnames(tp_fp, tp_fp_header) 
  colnames(tp_fp) <- make.unique(colnames(tp_fp))
  tp_fp <- tp_fp %>%
    filter(cnv!="neutral"&cnv.1!="neutral") %>%
    mutate(result_TF = ifelse(cnv==cnv.1, "TP", "FP")) %>%
    mutate(result_BP_start = ifelse(result_TF=="TP"& abs(start-start.1) < 100, "TP", "FP" )) %>%
    mutate(result_BP_end = ifelse(result_TF=="TP"& abs(end-end.1) < 100, "TP", "FP" ))
  if (!exists("temp_fp") || !file.exists(temp_fp) || file.size(temp_fp) == 0){
    fp <- data.table(matrix(ncol = length(colnames(tp_fp)), nrow = 0))
    fp <- setnames(fn, colnames(tp_fp))  # 设定列名
  } else {
    fp <- fread(temp_fp)
    fp <- setnames(fp, test_col)
    append_col_fp <- setdiff(colnames(tp_fp), test_col)
    for (col in append_col_fp) {
      set(fp, j = col, value = NA)
    }
    fp <- fp %>%
      filter(cnv != "neutral") %>%
      mutate(result_TF = "FP",
             result_BP_start = "FP",
             result_BP_end = "FP")
  }
  if (!exists("temp_fn") || !file.exists(temp_fn) || file.size(temp_fn) == 0) {
    fn <- data.table(matrix(ncol = length(colnames(tp_fp)), nrow = 0))
    fn <- setnames(fn, colnames(tp_fp))
  } else {
    fn <- fread(temp_fn)
    fn <- setnames(fn, true_col)
    append_col_fp <- setdiff(colnames(tp_fp), true_col)
    for (col in append_col_fp) {
      set(fn, j = col, value = NA)
    }
    fn <- fn %>%
      filter(cnv != "neutral") %>%
      mutate(result_TF = "FN",
             result_BP_start = "FN",
             result_BP_end = "FN")
  }
  
  merged_dt <- rbindlist(list(tp_fp, fn, fp), use.names = TRUE, fill = TRUE) %>%
    select(-result_BP_start, -result_BP_end)
  dt <- as.data.table(merged_dt)
  
  # 统计 TP, FP, FN 的数量
  tp_count <- nrow(dt[result_TF == "TP"])
  fp_count <- nrow(dt[result_TF == "FP"])
  fn_count <- nrow(dt[result_TF == "FN"]) 
  recall <- tp_count / (tp_count + fn_count)   # Sensitivity
  precision <- tp_count / (tp_count + fp_count) # Positive Predictive Value (PPV)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  metrics <- data.table(
    TP = tp_count,
    FP = fp_count,
    FN = fn_count,
    Recall = round(recall, 4),
    Precision = round(precision, 4),
    F1_Score = round(f1_score, 4)
  )
  if (!dir.exists(opt$out_dir)) {
    dir.create(opt$out_dir, recursive = TRUE)
  }
  
  # 保存结果
  out_dir <- opt$out_dir
  out_dir <- file.path(out_dir, "event")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)  # recursive=TRUE 确保父目录存在
  }
  merged_output <- file.path(out_dir, "event_results.tsv")
  metrics_output <- file.path(out_dir, "event_summary.tsv")
  
  fwrite(merged_dt, merged_output, sep = "\t", quote = FALSE, row.names = FALSE)
  fwrite(metrics, metrics_output, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Results saved to:\n")
  cat(paste0("- ", merged_output, "\n"))
  cat(paste0("- ", metrics_output, "\n"))


} else if (func == "breakpoint") {
  breakpioint <- standard %>%
    mutate(last_cnv = lag(cnv, default = "neutral")) %>%
    group_by(chr) %>%
    mutate(continuous = ifelse( start-lag(end, default = 1) < 1, 1, 0)) %>%
    ungroup() %>%
    mutate(left_cnv = ifelse(continuous == 1, last_cnv, "neutral")) %>%
    mutate(continuous_end = lead(continuous, default = 0)) %>%
    mutate(next_cnv = lead(cnv, default = "neutral")) %>%
    mutate(right_cnv = ifelse(continuous_end == 1, next_cnv, "neutral"))
  all_start <- breakpioint %>%
    select(left_cnv, chr, start, cnv)
  all_start <- setnames(all_start, c("left_cnv", "chr", "position", "right_cnv"))
  all_end <- breakpioint %>%
    select(cnv, chr, end, right_cnv)
  all_end <- setnames(all_end, c("left_cnv", "chr", "position", "right_cnv"))
  all <- rbind(all_start, all_end) %>%
    group_by(chr, position) %>%
    summarise(across(everything(), first)) %>%
    filter(left_cnv != right_cnv)
  
  break_true <- true_pos %>%
    mutate(last_cnv = lag(cnv, default = "neutral")) %>%
    group_by(chr) %>%
    mutate(continuous = ifelse( start-lag(end, default = 1) < 1, 1, 0)) %>%
    ungroup() %>%
    mutate(left_cnv = ifelse(continuous == 1, last_cnv, "neutral")) %>%
    mutate(continuous_end = lead(continuous, default = 0)) %>%
    mutate(next_cnv = lead(cnv, default = "neutral")) %>%
    mutate(right_cnv = ifelse(continuous_end == 1, next_cnv, "neutral"))
  all_start_true <- break_true %>%
    select(left_cnv, chr, start, cnv)
  all_start_true <- setnames(all_start_true, c("left_cnv", "chr", "position", "right_cnv"))
  all_end_true <- break_true %>%
    select(cnv, chr, end, right_cnv)
  all_end_true <- setnames(all_end_true, c("left_cnv", "chr", "position", "right_cnv"))
  all_true <- rbind(all_start_true, all_end_true) %>%
    group_by(chr, position) %>%
    summarise(across(everything(), first)) %>%
    filter(left_cnv != right_cnv)
  result <- fuzzy_full_join(
    all_true, all,
    by = c("chr" = "chr", "position" = "position"),  # 确保匹配的列名
    match_fun = list(`==`, function(x, y) abs(x - y) <= 50))  %>%
    mutate(result = case_when(
      !is.na(chr.x) & !is.na(chr.y) & left_cnv.x == left_cnv.y & right_cnv.x == right_cnv.y ~ "TP",
      is.na(chr.x) ~ "FP",
      is.na(chr.y) ~ "FN",
      (!is.na(chr.x) & !is.na(chr.y) & left_cnv.x != left_cnv.y) |
        (!is.na(chr.x) & !is.na(chr.y) & right_cnv.x != right_cnv.y) ~ "FP"
    ))
  dt_bp <- as.data.table(result)
  TP <- dt_bp[result == "TP", .N]  # 统计 TP 的数量
  FN <- dt_bp[result == "FN", .N]  # 统计 FN 的数量
  FP <- dt_bp[result == "FP", .N]  # 统计 FP 的数量（当前数据没有 FP，但代码仍保留计算）
  
  # 计算 Recall, Precision, F1-score
  recall <- TP / (TP + FN)
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))  # 避免除零错误
  F1 <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
  result_table <- data.table(
    TP = TP, FP = FP, FN = FN, 
    Recall = round(recall, 4), 
    Precision = round(precision, 4), 
    F1_Score = round(F1, 4)
  )
  # 输出结果
  cat("Recall:", recall, "\n")
  cat("Precision:", precision, "\n")
  cat("F1-score:", F1, "\n")
  out_dir <- file.path(out_dir, "breakpoint")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)  # recursive=TRUE 确保父目录存在
  }
  result_filename <- file.path(out_dir, "breakpoint_result.tsv")
  metrics_filename <- file.path(out_dir, "breakpoint_metrics.tsv")
  # 保存 result
  write_tsv(result, result_filename)
  
  # 保存 metrics 为 TSV
  write_tsv(result_table, metrics_filename)
  
  
  # 打印路径，确认保存成功
  cat("Result saved to:", result_filename, "\n")
  cat("Metrics saved to:", metrics_filename, "\n")
}

