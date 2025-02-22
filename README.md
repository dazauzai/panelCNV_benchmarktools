### README: Benchmarking CNV Detection Scripts

## Introduction
This script is designed for benchmarking copy number variation (CNV) detection methods. It provides functionality for standardizing data, evaluating CNV calls against true positive datasets, and generating benchmarking reports.

## Requirements
Before running the script, ensure you have the following dependencies installed:
- **Bash** (Linux/macOS environment)
- **R (>= 4.0.0)** with required packages
- **jq** (for JSON parsing)
- **Required R packages**: `data.table`, `dplyr`, `ggplot2`

## Usage

```bash
bash main.sh [OPTIONS]
```

### Options
| Option | Short | Description |
|--------|-------|-------------|
| `--seg` | `-s` | Path to the segment file |
| `--true` | `-T` | Path to the true positive CNV file |
| `--test` | `-t` | Path to the benchmark CNV file |
| `--tools` | `-l` | Path to tools directory |
| `--output` | `-o` | Output directory |
| `--name` | `-n` | Name for the output subdirectory |
| `--col` | `-c` | Column index for CN values (default: 4) |
| `--func` | `-f` | Function type: `gene`, `event`, `breakpoint`, `intersection_percent` |
| `--gender` | `-g` | Gender of the sample (XY/XX) |
| `--json` | `-j` | Path to a JSON file containing parameters |
| `--batch` | `-b` | Name of this batch |
| `--skip` | - | Skip standardization step |
| `--check` | - | Check environment dependencies |
| `--help` | `-h` | Show help message |

## Example Commands

### **1. Checking the environment**
```bash
bash main.sh --check
```

### **2. Running Standardization**
```bash
bash main.sh -t benchmark_data.txt -l tools_dir -o results -n test_run -f standarlization
```

### **3. Running Gene Benchmarking**
```bash
bash main.sh -t benchmark_data.txt -T true_cnv.txt -s segments.txt -o results -n test_run -f gene
```

### **4. Running Breakpoint Evaluation**
```bash
bash main.sh -t benchmark_data.txt -T true_cnv.txt -o results -n test_run -f breakpoint
```

### **5. Running Intersection Percent Analysis**
```bash
bash main.sh -t benchmark_data.txt -T true_cnv.txt -o results -n test_run -f intersection_percent
```

### **6. Running with a JSON Configuration**
Create a JSON file (`config.json`) with the following format:
```json
{
  "seg": "segments.txt",
  "true": "true_cnv.txt",
  "batch": "test_batch"
}
```
Then, run:
```bash
bash main.sh -j config.json -o results -n test_run -f gene
```

## Output Structure
The output directory will be structured as follows:
```
/output_directory/
    ├─ batch_name/
    │   ├─ run_name/
    │   │   ├─ gene_benchmark_results.txt
    │   │   ├─ breakpoint_analysis.txt
    │   │   ├─ intersection_percent_results.txt
    │   │   └─ config.json
```

## Troubleshooting
- **Error: "File not found"**
  - Ensure that the input files exist and that the paths are correctly specified.
- **Error: "Unknown function type"**
  - Ensure the `-f` option is set to one of: `gene`, `event`, `breakpoint`, `intersection_percent`.
- **Error: "jq command not found"**
  - Install `jq` using:
    ```bash
    sudo apt install jq  # Ubuntu/Debian
    brew install jq      # macOS
    ```

## Contact
For any issues, please reach out to **yiming@yourdomain.com** or submit a bug report in your project repository.

---

This README provides detailed instructions on how to use the benchmarking script efficiently. Let me know if you'd like additional clarifications or modifications! 🚀

