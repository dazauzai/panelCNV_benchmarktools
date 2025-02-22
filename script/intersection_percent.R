library(data.table)
library(dplyr)
library(optparse)
option_list <- list(
  make_option(c("-s", "--stand"), type = "character", default = NULL,
              help = "Path to the standard CNV file", metavar = "character"),
  make_option(c("-t", "--true"), type = "character", default = NULL,
              help = "Path to the true positive CNV file", metavar = "character"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "Output directory for results", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# 检查参数是否提供
if (is.null(opt$stand) || is.null(opt$true) || is.null(opt$out_dir)) {
  print_help(opt_parser)
  stop("Error: Missing required arguments. Please specify --stand, --true, and --out_dir.", call. = FALSE)
}
result <- opt$stand
true <- opt$true
out_dir <- opt$out_dir
result_dt <- fread(result, header = T)
true_pos <- fread(true, header = T)
true_pos <- setnames(true_pos, c("chr", "start", "end", "gene", "cnv", "cn"))
temp_inter <- tempfile(fileext = ".bed")
temp_fp <- tempfile(fileext = ".bed")
temp_fn <- tempfile(fileext = ".bed")
rm_header <- function(input_file, output_file = tempfile(fileext = ".bed")) {
  cmd <- sprintf("tail -n +2 %s > %s", input_file, output_file)
  system(cmd)
  return(output_file)
}
temp_rm_header_true <- rm_header(true)
temp_rm_header_test <- rm_header(result)
cmd_tp_fp <- sprintf("bedtools intersect -wo -a %s -b %s > %s", temp_rm_header_test, temp_rm_header_true, temp_inter)
cmd_fp <- sprintf("bedtools intersect -v -a %s -b %s > %s", temp_rm_header_test, temp_rm_header_true, temp_fp)
cmd_fn <- sprintf("bedtools intersect -v -b %s -a %s > %s", temp_rm_header_test, temp_rm_header_true, temp_fn)
system(cmd_tp_fp)
system(cmd_fp)
system(cmd_fn)
intersect <- temp_inter
intersect_fp <- temp_fp
intersect_fn <- temp_fn
header_result <- colnames(result_dt)
header_true <- colnames(true_pos)
interaction <- fread(intersect, header = F)
if (file.exists(intersect_fp) && file.info(intersect_fp)$size > 0) {
  interaction_fp <- fread(intersect_fp, header = FALSE)
  setnames(interaction_fp, header_result)
  interaction_fp <- interaction_fp %>%
    filter(cnv != "neutral") %>%
    mutate(length = end - start) %>%
    mutate(type = "FP")
} else {
  header_result_tp <- c(header_result, "length", "type")
  interaction_fp <- data.table(matrix(ncol = length(header_result_tp), nrow = 0))
  setnames(interaction_fp, header_result_tp)
}

# 读取 intersection_fn
if (file.exists(intersect_fn) && file.info(intersect_fn)$size > 0) {
  interaction_fn <- fread(intersect_fn, header = FALSE)
  setnames(interaction_fn, header_true)
  interaction_fn <- interaction_fn %>%
    filter(cnv != "neutral") %>%
    mutate(length = end - start) %>%
    mutate(type = "FN")
} else {
  header_true_tp <- c(header_true, "length", "type")
  interaction_fn <- data.table(matrix(ncol = length(header_true_tp), nrow = 0))
  setnames(interaction_fn, header_true_tp)
}
header_inter <- c(header_result, header_true, "intersect")
interaction <- setnames(interaction, header_inter)
colnames(interaction) <- make.unique(colnames(interaction))
interaction <- interaction %>%
  filter(cnv!="neutral", cnv.1 != "neutral") %>%
  mutate(length = ifelse(cnv==cnv.1, intersect, 0)) %>%
  mutate(type = ifelse(cnv==cnv.1, "TP", "FP"))
intersection <- as.data.table(interaction)
interaction_fp <- as.data.table(interaction_fp)
interaction_fn <- as.data.table(interaction_fn)
intersect_cols <- colnames(intersection)
bind <- list(interaction, interaction_fp, interaction_fn)
final_result <- rbindlist(bind, use.names=TRUE, fill=TRUE)
calculate_metrics <- function(final_result) {
  # 确保输入是 data.table
  final_result <- as.data.table(final_result)
  
  # 计算 TP, FP, FN 总长度
  tp <- final_result %>%
    filter(type == "TP") %>%
    summarise(total_length = sum(length, na.rm = TRUE)) %>%
    pull(total_length)
  
  fp <- final_result %>%
    filter(type == "FP") %>%
    summarise(total_length = sum(length, na.rm = TRUE)) %>%
    pull(total_length)
  
  fn <- final_result %>%
    filter(type == "FN") %>%
    summarise(total_length = sum(length, na.rm = TRUE)) %>%
    pull(total_length)
  
  # 计算 Recall, Precision, F1-score
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), NA)
  
  # 创建单行 data.table
  result_table_total <- data.table(
    FP = fp,
    TP = tp,
    FN = fn,
    Recall = recall,
    Precision = precision,
    F1_Score = f1_score
  )
  
  return(result_table_total)
}


final_result_gain <- final_result %>%
  filter((is.na(cnv.1) & cnv=="gain") | ( cnv.1 == "gain" & is.na(cnv)) | ( !is.na(cnv) & !is.na(cnv.1) & cnv==cnv.1 & cnv == "gain" ) | ( !is.na(cnv) & !is.na(cnv.1) & cnv!=cnv.1 & cnv.1 == "gain"))
final_result_loss <- final_result %>%
  filter((is.na(cnv.1) & cnv=="loss") | ( cnv.1 == "loss" & is.na(cnv)) | ( !is.na(cnv) & !is.na(cnv.1) & cnv==cnv.1 & cnv == "loss" ) | ( !is.na(cnv) & !is.na(cnv.1) & cnv!=cnv.1 & cnv.1 == "loss"))
result_table_total <- calculate_metrics(final_result)
result_table_gain <- calculate_metrics(final_result_gain)
result_table_loss <- calculate_metrics(final_result_loss)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 定义文件路径
out_dir <- file.path(out_dir, "intersect")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)  # recursive=TRUE 确保父目录存在
}
intersect_result_path <- file.path(out_dir, "intersect_result.bed")
intersect_result_gain_path <- file.path(out_dir, "intersect_result_gain.bed")
intersect_result_loss_path <- file.path(out_dir, "intersect_result_loss.bed")

intersect_matrix_total_path <- file.path(out_dir, "intersect_matrixs_total.bed")
intersect_matrix_gain_path <- file.path(out_dir, "intersect_matrixs_gain.bed")
intersect_matrix_loss_path <- file.path(out_dir, "intersect_matrixs_loss.bed")

# 写出数据
fwrite(final_result, intersect_result_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(final_result_gain, intersect_result_gain_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(final_result_loss, intersect_result_loss_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

fwrite(result_table_total, intersect_matrix_total_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(result_table_gain, intersect_matrix_gain_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(result_table_loss, intersect_matrix_loss_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 确认写出成功
message("✅ 所有文件已成功写入到 ", out_dir)

