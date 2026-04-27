# CARLIN序列分析工具使用手册

## 目录

1. [安装与配置](#安装与配置)
2. [快速开始](#快速开始)
3. [命令详解](#命令详解)
4. [输入文件格式](#输入文件格式)
5. [输出文件格式](#输出文件格式)
6. [参数调整指南](#参数调整指南)
7. [应用示例](#应用示例)
8. [故障排除](#故障排除)
9. [性能优化](#性能优化)

## 安装与配置

### 系统要求
- Python 3.7或更高版本
- 至少2GB可用内存（用于大型数据集）

### 安装步骤

1. **克隆或下载工具**
```bash
git clone <repository-url>
cd carlin-align-tool
```

2. **安装Python依赖**
```bash
pip install -r requirements.txt
```

3. **验证安装**
```bash
python carlin_tool.py --version
python carlin_tool.py --help
```

### 依赖说明
- **numpy** (>=1.19.0): 数值计算和矩阵操作
- **biopython** (>=1.79): FASTQ/FASTA文件解析

## 快速开始

### 示例1：基本工作流程

```bash
# 1. 将FASTQ转换为TSV格式
python carlin_tool.py convert fastq-to-tsv \
  --fastq example_data/test.fastq.gz \
  --output my_sample_queries.tsv \
  --sample-name my_sample

# 2. 批量序列比对
python carlin_tool.py align \
  --reference example_data/reference.fa \
  --queries my_sample_queries.tsv \
  --output my_sample_alignments.json \
  --format json
```

### 示例2：完整分析流程

```bash
#!/bin/bash
# 完整分析脚本示例

# 设置变量
FASTQ_FILE="data/reads.fastq.gz"
REFERENCE="data/reference.fasta"
SAMPLE_NAME="experiment1"
OUTPUT_PREFIX="results/${SAMPLE_NAME}"

# 创建输出目录
mkdir -p results

# 步骤1: FASTQ转换
echo "步骤1: 转换FASTQ文件..."
python carlin_tool.py convert fastq-to-tsv \
  --fastq "${FASTQ_FILE}" \
  --output "${OUTPUT_PREFIX}_queries.tsv" \
  --sample-name "${SAMPLE_NAME}"

# 步骤2: 序列比对
echo "步骤2: 批量序列比对..."
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${OUTPUT_PREFIX}_queries.tsv" \
  --output "${OUTPUT_PREFIX}_alignments.json" \
  --format json

# 步骤3: 生成统计摘要
echo "步骤3: 生成统计摘要..."
python -c "
import json
with open('${OUTPUT_PREFIX}_alignments.json') as f:
    data = json.load(f)

total = len(data)
successful = len([r for r in data if 'error' not in r])
failed = total - successful

print(f'总计序列: {total}')
print(f'成功比对: {successful}')
print(f'失败比对: {failed}')

if successful > 0:
    scores = [r['score'] for r in data if 'score' in r and r['score'] is not None]
    avg_score = sum(scores) / len(scores) if scores else 0
    print(f'平均比对得分: {avg_score:.2f}')
"

echo "分析完成！"
```

## 命令详解

### 主命令结构
```
carlin_tool.py <command> <subcommand> [选项]
```

### 转换命令 (convert)

#### fastq-to-tsv
将FASTQ文件转换为TSV格式。

**选项**：
- `--fastq`: 输入FASTQ文件（支持.gz压缩）
- `--output`: 输出TSV文件路径
- `--sample-name`: 样本名称（默认: "sample"）

**示例**：
```bash
python carlin_tool.py convert fastq-to-tsv \
  --fastq reads.fastq.gz \
  --output reads.tsv \
  --sample-name patient_01
```

#### fastq-to-fasta
将FASTQ文件转换为FASTA格式。

**选项**：
- `--fastq`: 输入FASTQ文件（支持.gz压缩）
- `--output`: 输出FASTA文件路径
- `--sample-name`: 样本名称（默认: "sample"）

**示例**：
```bash
python carlin_tool.py convert fastq-to-fasta \
  --fastq reads.fastq.gz \
  --output reads.fasta \
  --sample-name cell_line_A
```

### 比对命令 (align)

并行批量比对查询序列到参考序列。自动利用多CPU核心加速处理。

**必需选项**：
- `--reference`: 参考序列FASTA文件
- `--queries`: 查询序列文件（TSV或FASTA格式）
- `--output`: 输出文件路径

**比对参数**：
- `--match-score`: 匹配得分（默认: 2.0）
- `--mismatch-penalty`: 错配惩罚（默认: -3.0）
- `--gap-open`: gap开启惩罚（默认: -2.0）
- `--gap-extend`: gap延伸惩罚（默认: -0.1）
- `--global`: 使用全局比对（默认使用半全局比对）

**并行参数**：
- `--threads`, `-t`: 并行进程数（默认：自动使用所有CPU核心）

**输出选项**：
- `--format`: 输出格式，可选 `json`、`tsv` 或 `all`（同时输出两种格式，默认: json）
- `--report`: 生成突变分析报告，可选 `json` 或 `html`（默认: 不生成）
- `--report-output`: 报告输出路径（默认: 基于 `--output` 文件名自动生成）

**示例**：
```bash
# 使用默认参数，输出JSON
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json

# 输出TSV格式
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results.tsv \
  --format tsv

# 同时输出JSON和TSV
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results/prefix \
  --format all

# 比对并生成HTML格式分析报告
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --report html

# 比对并生成JSON格式报告，指定报告路径
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --report json \
  --report-output my_report.json

# 使用自定义参数和全局比对
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.fasta \
  --output results.tsv \
  --format tsv \
  --match-score 3.0 \
  --mismatch-penalty -5.0 \
  --gap-open -1.5 \
  --gap-extend -0.05 \
  --global
```

## 输入文件格式

### 1. 参考序列文件 (FASTA)
单条或多条序列的FASTA文件。如果包含多条序列，默认使用第一条。

**格式**：
```
>sequence_id
ACGTACGTACGTACGT...
```

**示例**：
```
>CARLIN_reference
TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACG...
```

### 2. 查询序列文件

#### TSV格式（推荐）
5列制表符分隔文件：

| 列名 | 说明 | 示例 |
|------|------|------|
| readName | 序列标识符 | sample1_seq1 |
| cellBC | 细胞条形码（虚拟） | sample1 |
| UMI | 唯一分子标识符（虚拟） | UMI1 |
| readCount | 该序列的观测次数 | 2 |
| seq | 序列字符串 | ACGTACGT... |

**示例文件**：
```
readName	cellBC	UMI	readCount	seq
sample1_seq1	sample1	UMI1	2	TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGC...
sample1_seq2	sample1	UMI2	1	TCTGTGTGGGAGGGCTAAGAGGCCGCCGGACTGC...
```

#### FASTA格式
支持元数据在头部中的FASTA文件。

**格式**：
```
>readName cellBC=sample1 UMI=UMI1 readCount=2
ACGTACGTACGT...
```

### 3. FASTQ文件
标准FASTQ格式，支持gzip压缩。

**格式**：
```
@read_id
ACGTACGT
+
IIIIIIII
```

## 输出文件格式

### JSON格式（默认）
完整的比对结果，包含所有详细信息。

**结构**：
```json
[
  {
    "readName": "sample1_seq1",
    "cellBC": "sample1",
    "UMI": "UMI1",
    "readCount": 2,
    "score": 659.0,
    "aligned_ref": "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGC...",
    "aligned_query": "TCTGTGTGGGAGGGCTAAGAGGCCGCCGGACTGC...",
    "stats": {
      "matches": 331,
      "mismatches": 1,
      "gaps_in_ref": 0,
      "gaps_in_query": 0,
      "gap_blocks_ref": [],
      "gap_blocks_query": [],
      "avg_gap_len_ref": 0,
      "avg_gap_len_query": 0,
      "alignment_length": 332,
      "similarity": 0.997,
      "identity": 0.997,
      "score": 659.0
    }
  }
]
```

### TSV格式
简化的表格格式，便于统计分析。

**列说明**：
- `readName`, `cellBC`, `UMI`, `readCount`: 输入序列信息
- `score`: 比对得分
- `matches`: 匹配碱基数
- `mismatches`: 错配碱基数
- `gaps_in_query`: 查询序列中的gap数
- `similarity`: 相似度（0-1）
- `aligned_ref`: 对齐后的参考序列（gap位置为`-`）
- `aligned_query`: 对齐后的查询序列（gap位置为`-`）
- `error`: 错误信息（如果比对失败）

> 使用 `--format all` 可同时输出JSON和TSV两种格式，无需重复运行。

### 分析报告 (--report)

使用 `--report json` 或 `--report html` 可生成突变分析报告，包含：

**报告内容**：
- **摘要统计**: 总序列数、总Reads数、成功/失败比对、编辑效率
- **突变统计**: 突变序列数、突变Reads数、未突变序列数
- **突变类型分布**（含序列数和Reads数）:
  - 仅点突变（替换）
  - 仅删除
  - 仅插入
  - 组合类型（插入+删除、插入+点突变等）
- **突变详细指标**:
  - 点突变总数
  - 插入/删除事件数
  - 平均插入/删除长度
  - 最大插入/删除长度
- **突变序列明细**: 每条突变序列的具体信息（最多100条）

**JSON报告示例**：
```json
{
  "tool": "CARLIN序列分析工具",
  "version": "2.1.0",
  "summary": {
    "total_sequences": 3284,
    "total_reads_all": 6055,
    "successful_alignments": 3276,
    "failed_alignments": 8,
    "total_reads_successful": 6041,
    "mutated_sequences": 1409,
    "unmutated_sequences": 1867,
    "mutated_reads": 2598,
    "editing_efficiency_pct": 43.01
  },
  "mutation_types": {
    "only_insertion": {"sequences": 12, "reads": 18},
    "only_deletion": {"sequences": 856, "reads": 1580},
    "only_substitution": {"sequences": 285, "reads": 524},
    "insertion_and_deletion": {"sequences": 97, "reads": 175},
    "insertion_and_substitution": {"sequences": 5, "reads": 9},
    "deletion_and_substitution": {"sequences": 148, "reads": 285},
    "insertion_deletion_substitution": {"sequences": 6, "reads": 7}
  },
  "mutation_stats": {
    "total_point_mutations": 542,
    "total_insertion_events": 124,
    "total_deletion_events": 1176,
    "avg_insertion_length": 3.42,
    "avg_deletion_length": 12.87,
    "max_insertion_length": 28,
    "max_deletion_length": 156
  }
}
```

**HTML报告**: 自包含的HTML文件，包含卡片式摘要、彩色表格和可视化布局，可直接在浏览器中打开查看。

## 参数调整指南

### 针对不同应用场景的参数建议

#### 1. CARLIN基因编辑分析
检测大片段缺失和点突变：
```bash
--match-score 2.0
--mismatch-penalty -3.0
--gap-open -1.5      # 降低开启惩罚，鼓励检测缺失
--gap-extend -0.05   # 降低延伸惩罚，鼓励长连续gap
```

#### 2. 严格序列验证
高严格度匹配，减少假阳性：
```bash
--match-score 3.0    # 提高匹配得分
--mismatch-penalty -5.0  # 增加错配惩罚
--gap-open -3.0      # 增加gap惩罚
--gap-extend -0.5    # 增加延伸惩罚
--global             # 使用全局比对
```

#### 3. 长序列比对
针对长序列（>500bp）优化：
```bash
--match-score 1.0    # 降低匹配得分权重
--gap-open -0.5      # 显著降低gap惩罚
--gap-extend -0.01   # 极低延伸惩罚
```

### 参数影响分析

| 参数 | 增加效果 | 减少效果 |
|------|----------|----------|
| match_score | 强化匹配，减少gap | 弱化匹配，增加gap |
| mismatch_penalty | 减少错配，增加gap | 增加错配，减少gap |
| gap_open | 减少gap，增加错配 | 增加gap，减少错配 |
| gap_extend | 减少长gap，增加短gap | 增加长gap，减少短gap |

## 应用示例

### 示例1：CARLIN编辑效率分析

```bash
#!/bin/bash
# CARLIN编辑效率分析脚本

REFERENCE="data/CARLIN_reference.fasta"
FASTQ="data/EPSC2_L1_1.fq.gz"
SAMPLE="EPSC2_L1"
OUTPUT_DIR="results"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 转换FASTQ
echo "转换FASTQ文件..."
python carlin_tool.py convert fastq-to-tsv \
  --fastq "${FASTQ}" \
  --output "${OUTPUT_DIR}/${SAMPLE}_queries.tsv" \
  --sample-name "${SAMPLE}"

# 批量比对（使用CARLIN优化参数）并生成HTML分析报告
echo "运行序列比对..."
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${OUTPUT_DIR}/${SAMPLE}_queries.tsv" \
  --output "${OUTPUT_DIR}/${SAMPLE}_alignments.json" \
  --gap-open -1.5 \
  --gap-extend -0.05 \
  --report html

echo "分析完成！查看报告: ${OUTPUT_DIR}/${SAMPLE}_alignments_report.html"
```

比对的JSON结果和HTML报告在同一目录下，报告中已包含编辑效率（约43%）和详细的突变统计信息。也可以单独使用 `--format all` 同时输出TSV格式用于自定义统计分析。
```

### 示例2：批量处理多个样本

```bash
#!/bin/bash
# 批量处理多个FASTQ文件

REFERENCE="reference.fasta"
SAMPLES=("sample1" "sample2" "sample3" "sample4")

for SAMPLE in \"${SAMPLES[@]}\"; do
    echo \"处理样本: ${SAMPLE}\"
    
    FASTQ=\"data/${SAMPLE}.fastq.gz\"
    OUTPUT_PREFIX=\"results/${SAMPLE}\"
    
    # 检查文件是否存在
    if [ ! -f \"${FASTQ}\" ]; then
        echo \"警告: 文件不存在 - ${FASTQ}\"
        continue
    fi
    
    # 转换和比对
    python carlin_tool.py convert fastq-to-tsv \
        --fastq \"${FASTQ}\" \
        --output \"${OUTPUT_PREFIX}_queries.tsv\" \
        --sample-name \"${SAMPLE}\"
    
    python carlin_tool.py align \
        --reference \"${REFERENCE}\" \
        --queries \"${OUTPUT_PREFIX}_queries.tsv\" \
        --output \"${OUTPUT_PREFIX}_alignments.json\"
    
    echo \"完成: ${SAMPLE}\"
done

echo \"所有样本处理完成！\"
```

## 故障排除

### 常见问题

#### 1. "ImportError: No module named 'Bio'"
**问题**：缺少BioPython库。
**解决**：
```bash
pip install biopython
```

#### 2. "文件未找到"错误
**问题**：文件路径不正确或文件不存在。
**解决**：
- 检查文件路径是否正确
- 使用绝对路径或相对路径
- 确保文件有读取权限

#### 3. "内存不足"错误
**问题**：处理大型数据集时内存不足。
**解决**：
- 分批处理数据
- 增加系统内存
- 使用更小的序列子集测试

#### 4. 比对结果异常
**问题**：参数设置不合理导致比对结果不符合预期。
**解决**：
- 调整比对参数（参考[参数调整指南](#参数调整指南)）
- 检查参考序列和查询序列格式
- 验证序列质量

#### 5. 处理速度慢
**问题**：大型数据集处理时间过长。
**解决**：
- 使用更高效的硬件
- 分批处理数据
- 考虑使用并行处理（见[性能优化](#性能优化)）

### 调试模式

对于复杂问题，可以添加调试输出：

```python
# 在carlin_tool.py中添加调试代码
import logging
logging.basicConfig(level=logging.DEBUG)
```

## 性能优化

### 1. 内置多进程并行加速（推荐）

本工具内置了多进程并行加速功能，无需额外配置即可利用所有CPU核心。

**基本用法**：
```bash
# 自动使用所有CPU核心
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json

# 手动指定进程数
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --threads 8
```

**并行加速效果**（以6055条序列为例）：
| 线程数 | 预计时间 | 加速比 |
|--------|----------|--------|
| 1（单线程） | 100%基准 | 1x |
| 2 | ~50% | ~2x |
| 4 | ~25% | ~4x |
| 8 | ~15% | ~6-7x |

> 注：加速比受CPU核心数和序列长度影响。序列越长，并行效率越高。

### 2. 分批处理大型数据集

```bash
#!/bin/bash
# 分批处理大型TSV文件

REFERENCE="reference.fasta"
INPUT_TSV="large_queries.tsv"
BATCH_SIZE=1000
OUTPUT_DIR="batch_results"

mkdir -p "${OUTPUT_DIR}"

# 分割大文件
split -l ${BATCH_SIZE} "${INPUT_TSV}" "${OUTPUT_DIR}/batch_"

# 分批处理
for BATCH in "${OUTPUT_DIR}"/batch_*; do
    BATCH_NAME=$(basename "${BATCH}")
    OUTPUT_FILE="${OUTPUT_DIR}/${BATCH_NAME}_alignments.json"
    
    echo "处理批次: ${BATCH_NAME}"
    
    python carlin_tool.py align \
        --reference "${REFERENCE}" \
        --queries "${BATCH}" \
        --output "${OUTPUT_FILE}"
done

# 合并结果
echo "合并结果..."
# 合并JSON文件的代码...
```

### 3. 使用外部并行工具

```bash
#!/bin/bash
# 使用GNU Parallel进行并行处理

REFERENCE="reference.fasta"
SAMPLES=("sample1" "sample2" "sample3" "sample4")

# 定义处理函数
process_sample() {
    SAMPLE=$1
    FASTQ="data/${SAMPLE}.fastq.gz"
    OUTPUT="results/${SAMPLE}_alignments.json"
    
    python carlin_tool.py convert fastq-to-tsv \
        --fastq "${FASTQ}" \
        --output "temp/${SAMPLE}.tsv" \
        --sample-name "${SAMPLE}"
    
    python carlin_tool.py align \
        --reference "${REFERENCE}" \
        --queries "temp/${SAMPLE}.tsv" \
        --output "${OUTPUT}"
    
    echo "完成: ${SAMPLE}"
}

export -f process_sample

# 创建临时目录
mkdir -p temp results

# 并行处理（最多4个同时进行）
parallel -j 4 process_sample ::: "${SAMPLES[@]}"

# 清理临时文件
rm -rf temp
```

### 4. 内存使用优化

对于极大序列：
- 使用`--stats-only`模式（如果实现）
- 分批处理序列
- 考虑使用磁盘缓存

## 高级功能

### 自定义统计计算

```python
#!/usr/bin/env python3
"""
自定义统计计算示例
"""

import json
import pandas as pd
import numpy as np

def analyze_alignments(json_file):
    """分析比对结果"""
    
    with open(json_file) as f:
        data = json.load(f)
    
    # 转换为DataFrame
    df = pd.DataFrame([
        {
            'readName': r['readName'],
            'score': r.get('score', np.nan),
            'matches': r['stats']['matches'] if r.get('stats') else np.nan,
            'mismatches': r['stats']['mismatches'] if r.get('stats') else np.nan,
            'gaps': r['stats']['gaps_in_query'] if r.get('stats') else np.nan,
            'similarity': r['stats']['similarity'] if r.get('stats') else np.nan,
            'has_error': 'error' in r
        }
        for r in data
    ])
    
    # 基本统计
    print(f"总序列数: {len(df)}")
    print(f"成功比对: {len(df[~df['has_error']])}")
    print(f"失败比对: {len(df[df['has_error']])}")
    
    if len(df[~df['has_error']]) > 0:
        valid_df = df[~df['has_error']]
        
        print(f"\n比对得分统计:")
        print(f"  平均得分: {valid_df['score'].mean():.2f}")
        print(f"  最高得分: {valid_df['score'].max():.2f}")
        print(f"  最低得分: {valid_df['score'].min():.2f}")
        
        print(f"\n突变统计:")
        mutated = valid_df[(valid_df['mismatches'] > 0) | (valid_df['gaps'] > 0)]
        print(f"  突变序列数: {len(mutated)}")
        print(f"  突变比例: {len(mutated)/len(valid_df)*100:.1f}%")
        
        print(f"\nGap分布:")
        gap_counts = valid_df['gaps'].value_counts().sort_index()
        for gaps, count in gap_counts.items():
            print(f"  {gaps}个gap: {count}序列 ({count/len(valid_df)*100:.1f}%)")
    
    return df

if __name__ == "__main__":
    results_df = analyze_alignments("alignments.json")
    results_df.to_csv("alignment_statistics.csv", index=False)
```

## 技术支持

### 获取帮助
```bash
# 查看所有命令
python carlin_tool.py --help

# 查看特定命令帮助
python carlin_tool.py convert --help
python carlin_tool.py align --help
```

### 报告问题
遇到问题时，请提供：
1. 使用的命令和参数
2. 错误信息全文
3. 输入文件示例（如果可能）
4. 系统环境信息

### 版本信息
```bash
python carlin_tool.py --version
```

---

*本使用手册最后更新：2026年4月23日（v2.1.0：新增 --format all、--report 分析报告）*