    #!/bin/bash

    # 获取脚本所在目录
    script_dir=$(dirname "$(realpath "$0")")
    script=${script_dir}/script
    temp=${script_dir}/temp
    c=4  # 默认列索引
    skip="F"
    check="F"
    # 帮助信息
    show_help() {
        echo "Usage: $0 [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  -s, --seg       Path to segment file"
        echo "  -T, --true      Path to true positive CNV file"
        echo "  -t, --test      Path to benchmark CNV file"
        echo "  -l, --tools     Path to tools directory"
        echo "  -o, --output    Output directory"
        echo "  -n, --name      Name for output subdirectory"
        echo "  -c, --col       Column index for CN values (default: 4)"
        echo "  -f, --func      Function type (e.g., 'gene' or 'json')"
        echo "  -g, --gender    Gender of the sample (XY/XX)"
        echo "  -j, --json      Path to JSON file containing parameters"
        echo "  -b, --batch     name of this batch"
        echo "  --skip  skip the standardisation"
        echo "--check check the env for this tools"
        echo "  -h, --help      Show this help message and exit"
        exit 0
    }

    # 解析命令行参数
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            -s|--seg) s="$2"; shift ;;
            -T|--true) T="$2"; shift ;;
            -t|--test) t="$2"; shift ;;
            -l|--tools) l="$2"; shift ;;
            -o|--output) o="$2"; shift ;;
            -n|--name) n="$2"; shift ;;
            -c|--col) c="$2"; shift ;;
            -f|--func) f="$2"; shift ;;
            -g|--gender) g="$2"; shift ;;
            -b| --batch) b="$2"; shift ;;
            -j|--json) json_file="$2"; shift ;;
            --skip) skip="T" ;;   # 如果提供 --skip，则设置 skip 为 "T"
            --check) check="T" ;;
            -h|--help) show_help ;;  # 显示帮助信息
            *) echo "Unknown parameter: $1"; show_help ;;
        esac
        shift
    done
    if [[ "$check" == "T" ]]; then
        Rscript ${script}/check_env.R
        exit
    fi
    # **如果 `-j` 提供，解析 JSON 并赋值**
    if [[ -n "$json_file" ]]; then
        if [[ ! -f "$json_file" ]]; then
            echo "Error: JSON file $json_file not found."
            exit 1
        fi

        # 从 JSON 文件读取参数
        s=$(jq -r '.seg' "$json_file")
        T=$(jq -r '.true' "$json_file")
        b=$(jq -r '.batch' "$json_file")

    fi

    # 检查必须的参数是否提供
    MISSING_PARAMS=()

    # **如果没有提供 JSON，则检查这些参数**
    if [[ -z "$json_file" ]]; then
        if [[ "$func" == "gene" ]]; then
            [[ -z "$s" ]] && MISSING_PARAMS+=("-s|--seg")
        fi 
        [[ -z "$T" ]] && MISSING_PARAMS+=("-T|--true")
    fi

    # **无论是否提供 JSON，这些参数都必须有**
    if [[ "$skip" == "F" ]]; then
        [[ -z "$l" ]] && MISSING_PARAMS+=("-l|--tools")
    fi
    [[ -z "$o" ]] && MISSING_PARAMS+=("-o|--output")
    [[ -z "$g" ]] && MISSING_PARAMS+=("-g|--gender")
    [[ -z "$f" ]] && MISSING_PARAMS+=("-f|--func")
    [[ -z "$t" ]] && MISSING_PARAMS+=("-t|--test")
    [[ -z "$f" ]] && MISSING_PARAMS+=("-f|--func")
    [[ -z "$n" ]] && MISSING_PARAMS+=("-n|--name")
    # 如果仍有缺失参数，报错并退出
    if [[ ${#MISSING_PARAMS[@]} -gt 0 ]]; then
        echo "Error: Missing required parameters: ${MISSING_PARAMS[*]}"
        show_help
    fi

    # **检查文件是否存在（只有 JSON 模式才检查）**
    if [[ -z "$json_file" && "$func" == "gene" ]]; then
        for FILE in "$s" "$T" "$t"; do
            if [[ ! -f "$FILE" ]]; then
                echo "Error: File $FILE not found."
                exit 1
            fi
        done
    fi
    Rscript ${script}/check_env.R
    # **确保输出目录存在**
    OUTDIR="$o/$b/$n"
    mkdir -p "$OUTDIR"
    func=${f}
    # **如果 `-f json`，创建 JSON 并退出**
    # **执行任务**
    if [[ "$skip" == "F" ]]; then
       echo "Rscript ${script}/standarlization.R -i $t -o ${temp}/${n} -f ${l} -g $g -c $c"
       Rscript ${script}/standarlization.R -i $t -o ${temp}/${n} -f ${l} -g $g -c $c       
    else
        cp $t ${temp}/${n}
    fi
    if [[ "$f" == "standarlization" || "$f" == "stand" || "$f" == "std" ]]; then
	exit
    fi
    stand_z=${temp}/${n}
    if [[ "$f" == "json" ]]; then
        # **创建 JSON 配置文件**
        json_output="$OUTDIR/config_${b}.json"
        jq -n \
            --arg seg "$s" \
            --arg true "$T" \
            --arg batch "$b" \
            '{seg: $seg, true: $true, batch: $batch}' > "$json_output"

        echo "✅ JSON config saved at: $json_output"
        exit 0
    elif [[ "$f" == "gene" ]]; then
        echo "Rscript ${script}/gene_benchmark.R -t ${stand_z} -T ${T} -s $s -o $OUTDIR"
        Rscript ${script}/gene_benchmark.R -t ${stand_z} -T ${T} -s $s -o $OUTDIR
    elif [[ "$f" == "event" ]]; then
        echo "Rscript ${script}/event.R -s ${stand_z} -t ${T} -o $OUTDIR -f "event""
        Rscript ${script}/event.R -s ${stand_z} -t ${T} -o $OUTDIR -f "event"
    elif [[ "$f" == "breakpoint" || "$f" == "BP" ]]; then
        echo "Rscript ${script}/event.R -s ${stand_z} -t ${T} -o $OUTDIR -f "breakpoint""
        Rscript ${script}/event.R -s ${stand_z} -t ${T} -o $OUTDIR -f "breakpoint"
    elif [[ "$f" == "intersection_percent" || "$f" == "intersection" || "$f" == "int" ]]; then
        echo "Rscript ${script}/intersection_percent.R -s ${stand_z} -t ${T} -o $OUTDIR"
        Rscript ${script}/intersection_percent.R -s ${stand_z} -t ${T} -o $OUTDIR
    else
        # **如果是未知的 func，则报错**
        echo "❌ Error: Unknown function type '$f'. Allowed values: 'json', 'gene'."
        exit 1
    fi


    echo "✅ All checks passed. Output directory created: $OUTDIR"
