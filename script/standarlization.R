library(data.table)
library(dplyr)
library(optparse)
standard_name = c("chr", "start", "end", "cn")
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file path", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file path", metavar = "character"),
  make_option(c("-c", "--cn_col"), type = "integer", default = 4,
              help = "Column index for CN values (default: 4)", metavar = "integer"),
  make_option(c("-f", "--func"), type = "character", default = "else",
              help = "Function type (default: 'else')", metavar = "character"),
  make_option(c("-g", "--gender"), type = "character", default = "XY",
              help = "gender of the sample(num of X)", metavar = "character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 赋值解析后的参数
input <- opt$input
output <- opt$output
cn_col <- opt$cn_col
func <- opt$func
gender <- opt$gender
if (is.null(opt$gender)) {
  warning("Warning: No gender (-g) specified. Defaulting to 'XY'.")
}
if (func == "cnvkit") {
# 读取文件
  cnvkit <- fread(input, header = T)
  
  # 定义一个函数来检测列是否存在
  check_column <- function(df, col_name) {
    if (!col_name %in% colnames(df)) {
      stop(paste("Error: Column", col_name, "not found in the dataset, please check the header and check if cns file has been processed by cnvkit.py call"))
    }
  }
  
  # 检测是否存在 'cn' 列
  check_column(cnvkit, "cn")
  
  # 选择指定列
  cnvkit <- cnvkit %>%
    select(chromosome, start, end, cn)
  
  standardized <- setnames(cnvkit, standard_name)
} else if (func == "cnv-z") {
  cnv_z <- fread(input, header = F) %>%
   select(V1, V2, V3, V10)
  standardized <- setnames(cnv_z, standard_name)
} else if (func == "control_freec") {
  control_freec <- fread(input, header = F) %>%
    select(V1, V2, V3, V4)
  standardized <- setnames(control_freec, standard_name)
} else if (func == "decon") {
} else if (func == "panel_cn"){
  panel_cn <- fread(input, header = T) %>%
    select(Chr, Start, End, CN) %>%
    mutate(CN=gsub("CN", "", CN))
  standardized <- setnames(panel_cn, standard_name)
} else {
  else_file <- fread(input, header = F)
  cols <- c("V1", "V2", "V3", paste0("V", as.numeric(c)))
  else_file <- else_file %>%
    select(all_of(cols))
  standardized <- setnames(else_file, standard_name)
}
if (func == "decon"){
 decon <-  fread(input, header = T) %>%
   group_by(CNV.ID) %>%
   summarise(across(everything(), first)) %>%
   ungroup() %>%
   select(Chromosome, Start, End, CNV.type) %>%
   mutate(CNV.type = ifelse(CNV.type=="deletion", "loss", "gain"))
 standardized <- setnames(decon, c("chr", "start", "end","cnv"))
 check_chr <- standardized$chr[1]
 if (nchar(check_chr)==1) {
   warning('Detected the column "chr" lacks prefix "chr", automatically adding')
   standardized <- standardized %>%
     mutate(chr=paste0("chr", chr))
 }
} else {
  check_chr <- standardized$chr[1]
  if (nchar(check_chr)==1) {
    warning('Detected the column "chr" lacks prefix "chr", automatically adding')
    standardized <- standardized %>%
      mutate(chr=paste0("chr", chr))
  }
  standardized <- standardized %>%
    mutate(cnv = case_when(
      # XY 男性，非 X/Y 染色体
      gender == "XY" & !(chr %in% c("chrX", "chrY")) & cn == 2 ~ "neutral",
      gender == "XY" & !(chr %in% c("chrX", "chrY")) & cn < 2 ~ "loss",
      gender == "XY" & !(chr %in% c("chrX", "chrY")) & cn > 2 ~ "gain",
      
      # XY 男性，X 或 Y 染色体
      gender == "XY" & chr %in% c("chrX", "chrY") & cn == 1 ~ "neutral",
      gender == "XY" & chr %in% c("chrX", "chrY") & cn < 1 ~ "loss",
      gender == "XY" & chr %in% c("chrX", "chrY") & cn > 1 ~ "gain",
      
      # XX 女性，任何染色体
      gender == "XX" & cn == 2 ~ "neutral",
      gender == "XX" & cn < 2 ~ "loss",
      gender == "XX" & cn > 2 ~ "gain",
      
      # 其他情况标记为 "ERROR"
      TRUE ~ "ERROR"
    ))
  
  # **检查是否有 "ERROR"，如果有就报错**
  if (any(standardized$cnv == "ERROR")) {
    stop("Error: Some rows could not be classified. Check 'gender', 'chr', and 'cn'.")
  }
}
write.table(standardized, file = output, col.names = T, row.names = F, sep = "\t", quote = F)
