library(data.table)
library(dplyr)
library(tidyr)
library(optparse)

standard_name = c("chr", "start", "end", "cn")

option_list <- list(
  make_option(c("-T", "--true"), type = "character", default = NULL,
              help = "Path to true positive CNV file", metavar = "character"),
  make_option(c("-t", "--test"), type = "character", default = NULL,
              help = "Path to benchmark CNV file", metavar = "character"),
  make_option(c("-s", "--seg"), type = "character", default = NULL,
              help = "Path to segment file with gene annotations", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 赋值解析后的参数
true <- opt$true
stand <- opt$test
seg <- opt$seg
output_dir <- opt$output_dir
true_pos <- fread(true, header = T)
standard <- fread(stand, header = T)

# 确保 standard 数据中有 "cnv" 列
if (!"cnv" %in% colnames(standard)) {
  stop("Error: no cnv col in test samples, please check your standardization")
}

segment <- fread(seg)
temp_intersect <- tempfile(fileext = ".bed")
temp_rmheader <- tempfile(fileext = ".bed")

cmd_rmheader <- sprintf("tail -n +2 %s > %s", stand, temp_rmheader)
cmd_intersect <- sprintf("bedtools intersect -wo -a %s -b %s > %s", temp_rmheader, seg, temp_intersect)
system(cmd_rmheader)
system(cmd_intersect)

intersect <- fread(temp_intersect)
segment <- setnames(segment, c('chr', "start", "end", "gene", "probes"))
segment_col <- colnames(segment)
stand_col <- colnames(standard)
intersect_col <- c(stand_col, segment_col, "intersect")
intersect <- setnames(intersect, intersect_col)
colnames(intersect) <- make.unique(colnames(intersect))

test_gene <- intersect %>%
  filter(cnv!="neutral") %>%
  select(gene, cnv, probes) %>%
  group_by(gene, cnv) %>%
  summarise(across(everything(), first)) %>%
  ungroup()

segment_gene <- segment %>%
  select(gene, probes)

true_gene <- true_pos %>%
  left_join(segment_gene, by="gene") %>%
  select(gene, cnv, probes) %>%
  group_by(gene, cnv) %>%
  summarise(across(everything(), first)) %>%
  ungroup() %>%
  full_join(test_gene, by = c("gene"))   %>%
  mutate(probes = coalesce(probes.x, probes.y)) %>%  # 选取不为 NA 的 probes
  select(-probes.x, -probes.y) %>%
  mutate(result = case_when(
    is.na(cnv.x) ~ "FP" ,
    is.na(cnv.y) ~ "FN" ,
    cnv.x == cnv.y ~ "TP" ,
    cnv.x != cnv.y ~ "FP"
  ))

summary_table <- true_gene %>%
  count(result) %>%
  pivot_wider(names_from = result, values_from = n, values_fill = 0)

# 确保 summary_table 中存在 TP, FP, FN，如果缺失则填充 0
summary_table <- summary_table %>%
  mutate(
    TP = ifelse("TP" %in% colnames(.), TP, 0),
    FP = ifelse("FP" %in% colnames(.), FP, 0),
    FN = ifelse("FN" %in% colnames(.), FN, 0)
  )

# 计算 Precision, Recall, F1-score，确保不会因 0 导致错误
summary_table <- summary_table %>%
  mutate(
    Precision = ifelse(TP == 0 & FP == 0, 0, TP / (TP + FP)),
    Recall = ifelse(TP == 0 & FN == 0, 0, TP / (TP + FN)),
    F1_Score = ifelse(TP == 0, 0, 
                      2 * (Precision * Recall) / (Precision + Recall))
  )

output_dir <- file.path(output_dir, "gene")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # recursive=TRUE 确保父目录存在
}

# 定义输出文件路径
true_gene_output <- file.path(output_dir, "true_gene_results.tsv")
summary_output <- file.path(output_dir, "summary_gene_results.tsv")

# 输出 true_gene
write.table(true_gene, file = true_gene_output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 输出 summary_table
write.table(summary_table, file = summary_output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 确保文件保存成功
message("Results saved to: ", true_gene_output)
message("Summary saved to: ", summary_output)
