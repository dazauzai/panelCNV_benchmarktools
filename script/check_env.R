# 指定需要的 R package
required_packages <- c("data.table", "dplyr", "optparse", "tidyr", "fuzzyjoin", "readr")

# 检查哪些 package 没有安装
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# 安装缺失的 package
if (length(missing_packages) > 0) {
  cat("以下 R 包未安装，正在安装: ", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org/")
} else {
  cat("✅ 所有 R 包都已安装，无需安装。\n")
}

# 加载所有包
lapply(required_packages, library, character.only = TRUE)
cat("✅ 所有 R 包已成功加载！\n")
