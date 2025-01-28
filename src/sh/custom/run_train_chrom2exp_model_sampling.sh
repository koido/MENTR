#!/bin/bash

set -eu

# デフォルト値の設定
IN_PATH="FANTOM5"
OUT="where you want to output files"
N_SAMPLING=1000
SEED_SAMPLING=12345
PY_SCRIPT_PATH="src/py3/custom/train_chrom2exp_model_sampling.py"
L1=0
L2=100
ETA=0.01
NUM_ROUND=100
BASE_SCORE=2
FILTER_STR="all"
PSEUDOCOUNT=0.0001
EARLY_STOPPING_ROUNDS=10
SEED=12345
THREADS=1
sample_ontology_File="resources/supp_table_10.sample_ontology_information.tsv.gz"
sample_ontology_ID="CL:0000034"

# ヘルプメッセージの表示
function help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  -i TEXT  入力パス [default: ${IN_PATH}]"
  echo "  -o TEXT  出力ディレクトリ [default: ${OUT}]"
  echo "  -n INT   サンプリング数 [default: ${N_SAMPLING}]"
  echo "  -S INT   サンプリングのシード値 [default: ${SEED_SAMPLING}]"
  echo "  -p TEXT  Pythonスクリプトのパス [default: ${PY_SCRIPT_PATH}]"
  echo "  -1 FLOAT L1正則化パラメータ [default: ${L1}]"
  echo "  -2 FLOAT L2正則化パラメータ [default: ${L2}]"
  echo "  -e FLOAT 学習率 [default: ${ETA}]"
  echo "  -r INT   ラウンド数 [default: ${NUM_ROUND}]"
  echo "  -b FLOAT ベーススコア [default: ${BASE_SCORE}]"
  echo "  -f TEXT  フィルター文字列 [default: ${FILTER_STR}]"
  echo "  -P FLOAT 疑似カウント [default: ${PSEUDOCOUNT}]"
  echo "  -E INT   アーリーストッピングラウンド [default: ${EARLY_STOPPING_ROUNDS}]"
  echo "  -s INT   一般的なシード値 [default: ${SEED}]"
  echo "  -t INT   スレッド数 [default: ${THREADS}]"
  echo "  -F TEXT  オントロジーファイル [default: ${sample_ontology_File}]"
  echo "  -I TEXT  オントロジーID [default: ${sample_ontology_ID}]"
  echo "  -h       このヘルプメッセージを表示"
  exit 1
}

# オプションの解析
while getopts "i:o:n:S:p:1:2:e:r:b:f:P:E:s:t:F:I:h" opt; do
  case $opt in
    i) IN_PATH="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    n) N_SAMPLING="$OPTARG" ;;
    S) SEED_SAMPLING="$OPTARG" ;;
    p) PY_SCRIPT_PATH="$OPTARG" ;;
    1) L1="$OPTARG" ;;
    2) L2="$OPTARG" ;;
    e) ETA="$OPTARG" ;;
    r) NUM_ROUND="$OPTARG" ;;
    b) BASE_SCORE="$OPTARG" ;;
    f) FILTER_STR="$OPTARG" ;;
    P) PSEUDOCOUNT="$OPTARG" ;;
    E) EARLY_STOPPING_ROUNDS="$OPTARG" ;;
    s) SEED="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    F) sample_ontology_File="$OPTARG" ;;
    I) sample_ontology_ID="$OPTARG" ;;
    h) help ;;
    \?) echo "Invalid option: -$OPTARG" >&2; help ;;
    :) echo "Option -$OPTARG requires an argument." >&2; help ;;
  esac
done

# パラメータの表示
echo "# == Parameters =="
echo "Input path: ${IN_PATH}"
echo "Output directory: ${OUT}"
echo "Number of sampling: ${N_SAMPLING}"
echo "Sampling seed: ${SEED_SAMPLING}"
echo "Python script path: ${PY_SCRIPT_PATH}"
echo "L1 regularization: ${L1}"
echo "L2 regularization: ${L2}"
echo "Learning rate: ${ETA}"
echo "Number of rounds: ${NUM_ROUND}"
echo "Base score: ${BASE_SCORE}"
echo "Filter string: ${FILTER_STR}"
echo "Pseudocount: ${PSEUDOCOUNT}"
echo "Early stopping rounds: ${EARLY_STOPPING_ROUNDS}"
echo "Seed: ${SEED}"
echo "Threads: ${THREADS}"
echo "Sample ontology file: ${sample_ontology_File}"
echo "Sample ontology ID: ${sample_ontology_ID}"
echo ""

# Pythonスクリプトの実行
python -u ${PY_SCRIPT_PATH} \
  --out_dir ${OUT} \
  --infile_dir ${IN_PATH} \
  --sample_ontology_ID ${sample_ontology_ID} \
  --sample_ontology_File ${sample_ontology_File} \
  --filterStr ${FILTER_STR} \
  --pseudocount ${PSEUDOCOUNT} \
  --num_round ${NUM_ROUND} \
  --early_stopping_rounds ${EARLY_STOPPING_ROUNDS} \
  --l1 ${L1} \
  --l2 ${L2} \
  --eta ${ETA} \
  --base_score ${BASE_SCORE} \
  --threads ${THREADS} \
  --seed ${SEED} \
  --n_sampling ${N_SAMPLING} \
  --seed_sampling ${SEED_SAMPLING}


