# CARLIN序列分析工具使用手册

## 目录

1. [安装与配置](#安装与配置)
2. [快速开始](#快速开始)
3. [命令详解](#命令详解)
4. [谱系示踪比对模式](#谱系示踪比对模式)
5. [输入文件格式](#输入文件格式)
6. [输出文件格式](#输出文件格式)
7. [参数调整指南](#参数调整指南)
8. [应用示例](#应用示例)
9. [故障排除](#故障排除)
10. [性能优化](#性能优化)

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

### 示例1：标准工作流程

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

### 示例2：谱系示踪分析

```bash
# 谱系示踪模式（自动推断cutsite位置，过滤假阳性突变）
python carlin_tool.py align \
  --reference example_data/reference.fa \
  --queries my_sample_queries.tsv \
  --output my_sample_lt_results.json \
  --lineage \
  --report html

# 参数微调
python carlin_tool.py align \
  --reference example_data/reference.fa \
  --queries my_sample_queries.tsv \
  --output my_sample_lt_results.json \
  --lineage \
  --cutsite-scale 1.0 \
  --far-scale 8.0 \
  --mutation-window 3
```

### 示例3：完整分析流程

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

# 步骤2: 谱系示踪比对并生成报告
echo "步骤2: 谱系示踪比对..."
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${OUTPUT_PREFIX}_queries.tsv" \
  --output "${OUTPUT_PREFIX}_alignments.json" \
  --lineage \
  --report html

echo "分析完成！查看报告: ${OUTPUT_PREFIX}_alignments_report.html"
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

并行批量比对查询序列到参考序列。支持标准比对和谱系示踪比对两种模式。

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

**谱系示踪参数**：
- `--lineage`: 启用谱系示踪比对模式
- `--cutsite-scale`: cutsite区域gap惩罚倍率（默认: 1.0，越小越容易开gap）
- `--flank-scale`: cutsite侧翼区gap惩罚倍率（默认: 2.0）
- `--far-scale`: 远离cutsite区域gap惩罚倍率（默认: 6.0，越大gap越难开启）
- `--flank-width`: cutsite侧翼范围bp（默认: 3）
- `--mutation-window`: 保留点突变的cutsite窗口半径bp（默认: 3）
- `--density-threshold`: mismatch密度阈值，超此阈值视为insertion（默认: 0.34）
- `--cutsites`: cutsite位置配置文件（JSON格式），不指定则自动推断标准CARLIN结构

**DP原生特征参数**（谱系模式推荐开启）：
- `--gap-exit-bonus`: gap→match转换额外惩罚（≤0，默认0.0，推荐-1.0），抑制indel碎片化
- `--short-match-window`: 短匹配区域阈值bp（默认0=关闭，推荐3-5），低于此长度的连续match打折
- `--short-match-discount`: 短匹配区域match_score折扣系数（0~1，默认1.0=不打折，推荐0.5）
- `--dense-mismatch-window`: 密集错配检测窗口bp（默认6）
- `--dense-mismatch-penalty`: 密集错配区域额外惩罚（≤0，默认0=关闭，推荐-2.0）
- `--homology-window`: 同源区域检测窗口bp（默认8）
- `--homology-penalty`: 同源区域match_score惩罚（≤0，默认0=关闭，推荐-1.0）
- `--isolated-base-penalty`: 孤立碱基匹配额外惩罚（≤0，默认0=关闭，推荐-2.0）

**矫正管线控制**：
- `--repeat-correction-mode`: 重复序列矫正模式，可选 `auto`（动态检测）、`hardcoded`（硬编码列表）、`off`（关闭），默认 `auto`
- `--disable-target-misalignment`: 关闭小片段跨靶点矫正
- `--disable-isolated-match-removal`: 关闭孤立匹配清除
- `--disable-dense-mismatch-correction`: 关闭后处理密集错配矫正
- `--disable-point-mutation-filtering`: 关闭点突变过滤

**并行参数**：
- `--threads`, `-t`: 并行进程数（默认：自动使用所有CPU核心）

**输出选项**：
- `--format`: 输出格式，可选 `json`、`tsv` 或 `all`（默认: json）
- `--report`: 生成突变分析报告，可选 `json` 或 `html`（默认: 不生成）
- `--report-output`: 报告输出路径（默认: 基于 `--output` 文件名自动生成）

**示例**：
```bash
# 使用默认参数，输出JSON
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output alignments.json

# 谱系示踪模式
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output lt_results.json \
  --lineage

# 谱系示踪模式+自定义参数
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output lt_results.json \
  --lineage \
  --cutsite-scale 1.0 \
  --far-scale 8.0 \
  --flank-width 5 \
  --mutation-window 4

# 谱系示踪模式+手动指定cutsite
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output lt_results.json \
  --lineage \
  --cutsites my_cutsites.json

# 同时输出JSON和TSV，生成HTML报告
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results/prefix \
  --format all \
  --report html
```

## 谱系示踪比对模式

### 适用场景

谱系示踪比对模式（`--lineage`）专为**多靶点基因编辑谱系示踪实验**设计，如CARLIN系统。其核心假设是：

1. 扩增子由**多个串联的target**组成，每个target包含**保守区 + cutsite**
2. 基因编辑（切割）主要发生在cutsite区域
3. 远离cutsite的mismatch很可能是**测序/PCR错误**而非真实突变
4. 连续高密度的mismatch实际上是**insertion事件**

### 工作原理

```
输入序列
    │
    ▼
┌─────────────────────────────────────┐
│ Step 1: 位置感知DP比对              │
│  ┌─────────────────────────────┐   │
│  │ cutsite处:    gap_open=-2.0 │   │
│  │ 侧翼±3bp:    gap_open=-4.0 │   │
│  │ 保守区域:    gap_open=-12.0│   │
│  └─────────────────────────────┘   │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 2: 高密度mismatch→indel转换    │
│  滑动窗口检测 >34% mismatch+≥2错配  │
│  → 将ref碱基替换为gap（insertion）  │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 3: 区域感知点突变过滤          │
│  cutsite±3bp内 → 保留               │
│  紧邻gap的突变 → 保留（例外规则）   │
│  其余 → 矫正为ref碱基               │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 4: 跨靶点重复序列矫正          │
│  检测11组重复序列（15bp~4bp）       │
│  将错误匹配到远端副本的碱基搬回     │
│  - ACTGCACGACAGTCG (T1↔T9, 15bp)   │
│  - ACTCGCG (T2↔T7, 7bp)            │
│  - ACAGTCG (T1↔T3↔T9, 7bp)         │
│  - GAGCGC / GCGACT / GATACG / ...   │
│  - GACGA (T1↔T3, 5bp)              │
│  - ACTA (T9↔T3, 4bp)               │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 5: 小片段跨靶点矫正            │
│  - TAGTAT → Target8                 │
│  - 单碱基A: Target9→Target1         │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│ Step 6: 孤立匹配清除                │
│  删除区域中的孤立单碱基匹配         │
│  合并连续deletion片段               │
└─────────────────────────────────────┘
    │
    ▼
输出结果
```

### cutsite位置配置

#### 自动推断（默认）
标准CARLIN扩增子（332bp）自动识别10个cutsite：
```bash
python carlin_tool.py align --lineage ...
```

#### 手动指定（JSON格式）
```json
{
  "cutsites": [
    {"name": "Target1", "start": 41, "end": 47},
    {"name": "Target2", "start": 68, "end": 74},
    {"name": "Target3", "start": 95, "end": 101},
    {"name": "Target4", "start": 122, "end": 128},
    {"name": "Target5", "start": 149, "end": 155},
    {"name": "Target6", "start": 176, "end": 182},
    {"name": "Target7", "start": 203, "end": 209},
    {"name": "Target8", "start": 230, "end": 236},
    {"name": "Target9", "start": 257, "end": 263},
    {"name": "Target10", "start": 284, "end": 290}
  ]
}
```
cutsite坐标为0-based、inclusive。

### 比对结果解读

谱系示踪模式在标准比对统计基础上增加以下指标：

| 指标 | 说明 |
|------|------|
| `n_mutations_corrected` | 被矫正的假阳性点突变数量 |
| `dense_regions_converted` | 是否触发了高密度mismatch→indel转换 |

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
完整的比对结果，包含所有详细信息。谱系示踪模式在后处理矫正后输出最终对齐结果。

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
      "matches": 329,
      "mismatches": 0,
      "gaps_in_ref": 3,
      "gaps_in_query": 0,
      "gap_blocks_ref": [3],
      "gap_blocks_query": [],
      "avg_gap_len_ref": 3.0,
      "avg_gap_len_query": 0,
      "alignment_length": 335,
      "similarity": 0.982,
      "identity": 1.0,
      "score": 659.0,
      "n_mutations_corrected": 0,
      "dense_regions_converted": false
    }
  }
]
```

**字段说明**：
- `aligned_ref`: 经后处理矫正后的参考序列对齐（含gap符号`-`）
- `aligned_query`: 经后处理矫正后的查询序列对齐（含gap符号`-`）
- `stats.n_mutations_corrected`: 被区域感知点突变过滤矫正的碱基数
- `stats.dense_regions_converted`: 是否触发了高密度mismatch→indel转换

**后处理矫正**：alin后的序列会依次经过跨靶点重复序列矫正、小片段跨靶点矫正和孤立匹配清除三个步骤，确保多靶点谱系示踪数据中因重复序列导致的错误匹配被纠正。

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

### 分析报告 (--report)

使用 `--report json` 或 `--report html` 可生成突变分析报告。

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

**示例**：
```bash
# 标准比对 + HTML报告
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results.json \
  --report html

# 谱系示踪比对 + HTML报告
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results.json \
  --lineage \
  --report html
```

## 参数调整指南

### 标准比对参数

| 参数 | 增加效果 | 减少效果 |
|------|----------|----------|
| match_score | 强化匹配，减少gap | 弱化匹配，增加gap |
| mismatch_penalty | 减少错配，增加gap | 增加错配，减少gap |
| gap_open | 减少gap，增加错配 | 增加gap，减少错配 |
| gap_extend | 减少长gap，增加短gap | 增加长gap，减少短gap |

### 谱系示踪参数

#### `--cutsite-scale`（默认：1.0）
cutsite内部的gap惩罚倍率。越小越容易在cutsite处开启gap：
- **0.5**: gap_open = -1.0（非常容易开gap）
- **1.0**: gap_open = -2.0（基准）
- **2.0**: gap_open = -4.0（保守，减少cutsite处gap）

#### `--far-scale`（默认：6.0）
保守/骨架区域的gap惩罚倍率。越大越抑制gap：
- **4.0**: gap_open = -8.0（相对宽松）
- **6.0**: gap_open = -12.0（基准，强烈抑制）
- **10.0**: gap_open = -20.0（极端抑制）

#### `--mutation-window`（默认：3）
保留点突变的窗口大小（bp）。cutsite两侧各保留此范围内的点突变：
- **2**: 更严格的突变过滤（cutsite±2bp）
- **3**: 默认（cutsite±3bp）
- **5**: 更宽松（cutsite±5bp）

#### `--density-threshold`（默认：0.34）
mismatch密度阈值。连续区域中mismatch比例超过此值则视为insertion：
- **0.25**: 更敏感，容易触发转换
- **0.34**: 默认平衡值（略高于1/3）
- **0.50**: 更保守，仅极高密度才转换

#### `--flank-scale`（默认：2.0）
cutsite侧翼区域（±flank_width范围内）的gap惩罚倍率：
- **1.0**: 与cutsite相同（宽松）
- **2.0**: 中等抑制（默认）
- **4.0**: 强抑制

### 不同应用场景的参数建议

#### 1. CARLIN谱系示踪分析（推荐）
```bash
--lineage \
--cutsite-scale 1.0 \
--flank-scale 2.0 \
--far-scale 6.0 \
--mutation-window 3 \
--density-threshold 0.34
```

#### 2. 严格突变检测（降低假阳性）
```bash
--lineage \
--cutsite-scale 1.0 \
--far-scale 8.0 \
--mutation-window 2 \
--density-threshold 0.40
```

#### 3. 宽松检测（捕获更多事件）
```bash
--lineage \
--cutsite-scale 0.8 \
--far-scale 4.0 \
--mutation-window 4 \
--density-threshold 0.30
```

#### 4. 标准比对（常规用途）
```bash
--match-score 2.0 \
--mismatch-penalty -3.0 \
--gap-open -2.0 \
--gap-extend -0.1
```

## 应用示例

### 示例1：谱系示踪编辑效率分析

```bash
#!/bin/bash
# 谱系示踪编辑效率分析脚本

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

# 谱系示踪比对 + HTML分析报告
echo "运行谱系示踪比对..."
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${OUTPUT_DIR}/${SAMPLE}_queries.tsv" \
  --output "${OUTPUT_DIR}/${SAMPLE}_results" \
  --format all \
  --lineage \
  --report html

echo "分析完成！查看报告: ${OUTPUT_DIR}/${SAMPLE}_results_report.html"
```

### 示例2：批量处理多个样本

```bash
#!/bin/bash
# 批量处理多个FASTQ文件

REFERENCE="reference.fasta"
SAMPLES=("sample1" "sample2" "sample3" "sample4")

for SAMPLE in "${SAMPLES[@]}"; do
    echo "处理样本: ${SAMPLE}"
    
    FASTQ="data/${SAMPLE}.fastq.gz"
    OUTPUT_PREFIX="results/${SAMPLE}"
    
    # 检查文件是否存在
    if [ ! -f "${FASTQ}" ]; then
        echo "警告: 文件不存在 - ${FASTQ}"
        continue
    fi
    
    # 转换和谱系示踪比对
    python carlin_tool.py convert fastq-to-tsv \
        --fastq "${FASTQ}" \
        --output "${OUTPUT_PREFIX}_queries.tsv" \
        --sample-name "${SAMPLE}"
    
    python carlin_tool.py align \
        --reference "${REFERENCE}" \
        --queries "${OUTPUT_PREFIX}_queries.tsv" \
        --output "${OUTPUT_PREFIX}_results.json" \
        --lineage
    
    echo "完成: ${SAMPLE}"
done

echo "所有样本处理完成！"
```

### 示例3：不同参数对比分析

```bash
#!/bin/bash
# 对比不同参数设置的效果

REFERENCE="reference.fasta"
QUERIES="queries.tsv"
OUTPUT_DIR="comparison"

mkdir -p "${OUTPUT_DIR}"

# 标准比对
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${QUERIES}" \
  --output "${OUTPUT_DIR}/standard.json"
echo "标准比对完成"

# 谱系示踪 - 默认参数
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${QUERIES}" \
  --output "${OUTPUT_DIR}/lineage_default.json" \
  --lineage
echo "谱系示踪(默认)完成"

# 谱系示踪 - 严格模式
python carlin_tool.py align \
  --reference "${REFERENCE}" \
  --queries "${QUERIES}" \
  --output "${OUTPUT_DIR}/lineage_strict.json" \
  --lineage \
  --far-scale 10.0 \
  --mutation-window 2
echo "谱系示踪(严格)完成"
```

### 示例4：自定义cutsite配置文件

```bash
# 创建cutsite配置文件
cat > my_cutsites.json << 'EOF'
{
  "cutsites": [
    {"name": "Target1", "start": 41, "end": 47},
    {"name": "Target2", "start": 68, "end": 74},
    {"name": "Target3", "start": 95, "end": 101},
    {"name": "Target4", "start": 122, "end": 128},
    {"name": "Target5", "start": 149, "end": 155},
    {"name": "Target6", "start": 176, "end": 182},
    {"name": "Target7", "start": 203, "end": 209},
    {"name": "Target8", "start": 230, "end": 236},
    {"name": "Target9", "start": 257, "end": 263},
    {"name": "Target10", "start": 284, "end": 290}
  ]
}
EOF

# 使用配置文件
python carlin_tool.py align \
  --reference ref.fasta \
  --queries queries.tsv \
  --output results.json \
  --lineage \
  --cutsites my_cutsites.json
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
- 谱系示踪模式下检查cutsite位置是否正确

#### 5. 谱系示踪模式下cutsite检测失败
**问题**：`--lineage` 无法自动推断cutsite位置。
**解决**：
- 确认参考序列为标准332bp CARLIN扩增子
- 使用 `--cutsites` 手动指定cutsite位置
- 检查参考序列是否正确（以Primer5开始，Primer3结束）

#### 6. 处理速度慢
**问题**：大型数据集处理时间过长。
**解决**：
- 使用更高效的硬件
- 分批处理数据
- 谱系示踪模式比标准模式慢约10-15%（额外后处理）

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

> 注：加速比受CPU核心数和序列长度影响。系统默认上限**12线程**以防止CPU过载和系统崩溃。

### 3. 向量化DP加速

谱系示踪算法v3.0引入了NumPy向量化DP递推，核心优化包括：
- **Iy和M按行向量化**：每行内所有列的Iy和M状态一次性计算
- **cumsum对角线密度计算**：使用NumPy cumsum沿矩阵对角线计算密度，替代O(m×n×window)三重循环
- **评分矩阵预计算**：在DP之前一次性计算完整的m×n评分矩阵

**加速效果**（332bp CARLIN，全部特性开启）：

| 优化项 | 优化前 | 优化后 | 加速比 |
|--------|--------|--------|--------|
| 单序列比对 | 1.60s | 0.37s | 4.3x |
| 批量500条 (12线程) | 34s | 19s | 1.8x |
| 全量23,430条 (12线程) | ~38min | ~15min | 2.5x |

### 4. 分批处理大型数据集

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
        --output "${OUTPUT_FILE}" \
        --lineage
done

echo "批量处理完成！"
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
        --output "${OUTPUT}" \
        --lineage
    
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

*本使用手册最后更新：2026年5月15日（v3.0.0：向量化DP加速 + DP原生特征 + 矫正管线控制）*|
