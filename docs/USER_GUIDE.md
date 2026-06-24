# User guide / 使用手册

## Installation / 安装

或从源码安装：

```bash
git clone https://github.com/sxgou/crisviper
cd crisviper
pip install -e .
```

或用 Docker：

```bash
docker build -t crisviper .
docker run --rm -v $(pwd)/data:/data crisviper align --reference /data/ref.fa --queries /data/q.tsv --output /data/out.json
```

验证安装：

```bash
crisviper --version
crisviper --help
```

---

## Quick start / 快速上手

### 1. 标准扩增子（单靶点编辑）

```bash
# FASTQ → TSV
crisviper convert fastq-to-tsv --fastq reads.fastq.gz --output reads.tsv --sample-name my_sample

# 比对
crisviper align --reference ref.fasta --queries reads.tsv --output results.json
```

输出 JSON：包含比对后的序列、得分和基本统计（匹配数、错配数、gap 数）。

### 2. 多靶点谱系示踪（如 CARLIN，10 targets）

```bash
crisviper align \
  --reference carlin_ref.fa \
  --queries reads.tsv \
  --output lineage_results \
  --format all \
  --lineage \
  --report html
```

`--lineage` 开启结构感知的梯度 gap 惩罚、自动检测 cutsite 位置。HTML 报告包含突变类型分布、长度分布、摘要统计和等位基因热图。

输出目录下还会生成 6 个 TSV 总结表格：`allele_frequency.tsv`、`per_target_editing.tsv`、`filter_reason.tsv`、`deletion_length.tsv`、`insertion_length.tsv`、`event_level_details.tsv`。

### 3. 单细胞 RNA-seq 谱系数据

```bash
crisviper convert fastq-to-tsv --fastq sc_reads.fastq.gz --output reads.tsv --sample-name liver_sample
crisviper align --reference ref.fa --queries reads.tsv --output results.json --lineage --report html
```

TSV 格式保留了 cellBC 和 UMI 列，下游的细胞级别分析可以直接使用比对结果。

### 4. 使用 DP 原生特征（谱系模式推荐）

```bash
crisviper align \
  --reference ref.fa \
  --queries reads.tsv \
  --output results.json \
  --lineage \
  --report html \
  --gap-exit-strength -1.0 \
  --short-match-window 3 --short-match-discount 0.5 \
  --dense-mismatch-penalty -2.0 \
  --homology-penalty -1.0 \
  --isolated-base-penalty -2.0
```

---

## Input formats / 输入格式

### Reference (FASTA)

单条或多条序列的 FASTA。多条时默认使用第一条。

### Queries (TSV)

五列制表符分隔：

| Column | Description / 说明 | Example |
|--------|-------------|---------|
| readName | 序列标识符 | sample1_seq47 |
| cellBC | 细胞条形码 | sample1 |
| UMI | 唯一分子标识符 | UMI_001 |
| readCount | 观测次数 | 3 |
| seq | 碱基序列 | ACGTACGT... |

### Queries (FASTA)

FASTA 头部携带元数据：

```
>seq_id cellBC=sample1 UMI=UMI_001 readCount=3
ACGTACGT...
```

### Cutsite config (JSON)

非标准扩增子结构（不是标准 332 bp CARLIN）时，提供 cutsite 配置：

```json
{"cutsites": [
  {"name": "Target1", "start": 41, "end": 47}
]}
```

---

## Commands reference / 命令参考

### `crisviper convert fastq-to-tsv`

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--fastq` | — | 输入 FASTQ（单端，支持 .gz） |
| `--fastq1` | — | 双端 R1 FASTQ（需与 `--fastq2` 共用） |
| `--fastq2` | — | 双端 R2 FASTQ（需与 `--fastq1` 共用） |
| `--output` | — | 输出 TSV 路径 |
| `--sample-name` | `sample` | 样本名称 |
| `--min-reads` | 1 | 最小 read 数过滤 |
| `--min-overlap` | 10 | 双端合并最小重叠长度 (bp) |
| `--max-mismatch-rate` | 20 | 重叠区域最大错配率 (%) |
| `--max-mismatch-diff` | 5 | 重叠区域最大绝对错配数 |
| `--require-qual` | 15 | 合并最低碱基质量 (Phred) |

### `crisviper convert fastq-to-fasta`

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--fastq` | — | 输入 FASTQ（单端，支持 .gz） |
| `--fastq1` | — | 双端 R1 FASTQ（需与 `--fastq2` 共用） |
| `--fastq2` | — | 双端 R2 FASTQ（需与 `--fastq1` 共用） |
| `--output` | — | 输出 FASTA 路径 |
| `--sample-name` | `sample` | 样本名称 |
| `--min-overlap` | 10 | 双端合并最小重叠长度 (bp) |
| `--max-mismatch-rate` | 20 | 重叠区域最大错配率 (%) |
| `--max-mismatch-diff` | 5 | 重叠区域最大绝对错配数 |
| `--require-qual` | 15 | 合并最低碱基质量 (Phred) |

### `crisviper align`

**基本选项：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--reference` | — | 参考序列 FASTA |
| `--queries` | — | 查询序列文件（TSV、FASTA 或 FASTQ） |
| `--fastq1` | — | 双端 R1 FASTQ 文件（需与 `--fastq2` 共用） |
| `--fastq2` | — | 双端 R2 FASTQ 文件（需与 `--fastq1` 共用） |
| `--output` | — | 输出路径（`--format all` 时作为前缀） |
| `--format` | `json` | 输出格式：`json`、`tsv`、`all` |
| `--config` | — | YAML 配置文件（靶标/扩增子结构和管道参数） |
| `--sample-name` | `sample` | FASTQ 输入时的样本标记名 |
| `--threads` / `-t` | 1 | 并行进程数 |
| `--chunk-size` | 500 | 每批处理的序列数；传 0 启用自动计算 |

**打分参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--match-score` | 2.0 | 匹配得分 |
| `--mismatch-penalty` | -3.0 | 错配惩罚 |
| `--gap-open` | -2.0 | Gap 开启惩罚 |
| `--gap-extend` | -0.1 | Gap 延伸惩罚 |

**谱系模式参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--lineage` | off | 开启谱系示踪比对（结构感知梯度 gap 惩罚） |
| `--min-scale` | 1.0 | 切割点处最低惩罚倍率（越低越易开 gap） |
| `--max-scale` | 6.0 | 保守区最高惩罚倍率 |
| `--cutsite-edge-scale` | 2.0 | Cutsite 边界惩罚倍率 |
| `--gradient-radius` | auto | 梯度有效半径 (bp)，省略则自动计算 |
| `--sub-window` | 3 | Cutsite 邻近保留窗口 (bp)，控制背景矫正和突变标注 |
| `--mismatch-density-threshold` | 0.34 | 密集错配检测密度阈值，超过则转为 indel |
| `--cutsites` | auto | Cutsite 配置文件路径（JSON 格式，省略则自动检测） |

**DP 原生特征（谱系模式推荐）：**

| Option | Default | Recommended | Description / 说明 |
|--------|---------|-------------|-------------|
| `--gap-exit-strength` | 0.0 | -1.0 | Gap→match 转换惩罚；合并碎片化 indel（≤0，0=关闭） |
| `--short-match-window` | 0 | 3 | 短匹配区域阈值 (bp)，0=关闭 |
| `--short-match-discount` | 1.0 | 0.5 | 短匹配区域得分折扣 (1.0=不打折) |
| `--dense-mismatch-window` | 6 | 6 | 密集错配检测窗口 (bp) |
| `--dense-mismatch-penalty` | 0.0 | -2.0 | 密集错配区域额外惩罚，0=关闭（≤0） |
| `--homology-window` | 8 | 8 | 同源区域检测窗口 (bp) |
| `--homology-penalty` | 0.0 | -1.0 | 同源区域惩罚，0=关闭（≤0） |
| `--isolated-base-penalty` | 0.0 | -2.0 | 孤立碱基匹配惩罚，0=关闭（≤0） |

**引物参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--primer5-len` | 23 | 5' 端引物长度 (bp) |
| `--primer3-len` | 33 | 3' 端引物长度 (bp) |
| `--primer5-threshold` | 19 | 5' 引物匹配碱基数阈值 |
| `--primer3-threshold` | 29 | 3' 引物匹配碱基数阈值 |

**配对端合并参数（仅 FASTQ 输入）：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--min-overlap` | 10 | 配对端合并最小重叠长度 (bp) |
| `--max-mismatch-rate` | 20 | 重叠区域最大错配率 (%) |
| `--max-mismatch-diff` | 5 | 重叠区域最大绝对错配数 |
| `--require-qual` | 15 | 合并最低碱基质量 (Phred) |

**Allele 过滤和输出参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--min-reads` | 1 | 输入侧最小 read 数阈值（预处理过滤用） |
| `--min-reads-sub` | 5 | 纯点突变 allele 最小 read 数阈值（inclusive，>=此值通过） |
| `--min-reads-indel` | 0 | 含 indel 的 allele 最小 read 数阈值（0=不过滤） |
| `--correct-bg-sub` | on | 启用背景点突变矫正 |
| `--no-correct-bg-sub` | — | 关闭背景点突变矫正 |
| `--keep-sub-indel-window` | 3 | 背景矫正时 indel 邻近保留窗口 (bp) |
| `--read-to-allele` | off | 输出 read→allele 映射表到输出文件夹（仅 FASTQ 输入） |

**报告选项：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--report` | — | 生成报告：`json` 或 `html` |
| `--report-output` | 从 `--output` 推断 | 报告路径（去除扩展名 + `_report` 后缀） |
| `--allele-top-n` | 50 | 报告中展示的 top N alleles |
| `--allele-window-start` | 0 | Allele 热图显示起始位置 |
| `--allele-window-end` | end of reference | Allele 热图显示结束位置（含，0-indexed） |

---

## Output / 输出格式

### JSON

比对结果数组。每个条目：

```json
{
  "readName": "sample1_seq47",
  "cellBC": "sample1",
  "UMI": "UMI_001",
  "readCount": 3,
  "score": 659.0,
  "aligned_ref": "TATGTGT--GGGAGGG...",
  "aligned_query": "TATGTGTCGGGAGGG...",
  "mutations": [
    {"type": "deletion", "ref_pos": 41, "ref_base": "A", "length": 2, "in_cutsite_window": true}
  ],
  "stats": {
    "matches": 329,
    "mismatches": 0,
    "gaps_in_ref": 2,
    "gaps_in_query": 0,
    "similarity": 0.982,
    "identity": 1.0,
    "has_mutation": true,
    "dense_regions_converted": false
  }
}
```

### TSV

简化的表格：readName, cellBC, UMI, readCount, score, matches, mismatches, gaps_in_query, similarity, aligned_ref, aligned_query, error。

### HTML report / 分析报告

`--report html` 生成的 HTML 报告包含：摘要统计、突变类型条形图、长度分布、top allele 表、等位基因热图和每条突变序列的突变标签。所有 JavaScript 内联，无外部依赖，离线可用。

### Summary tables / 总结表格

运行完成后（`--format all` 或默认 JSON 输出），在输出目录下自动生成 6 个 TSV 表格：

| 文件 | 内容 |
|------|------|
| `allele_frequency.tsv` | 按突变指纹聚合的 Allele 频率表：Rank, Allele, Mutation_Type, Sequences, Reads, Reads_Pct |
| `per_target_editing.tsv` | 每个 Target 的编辑类型统计：Total, Edited, Rate_Pct, Del, Ins, Sub, Avg_Mut_Length（使用 20bp 窗口） |
| `filter_reason.tsv` | 序列丢弃原因统计：Reason, Sequences, Reads |
| `deletion_length.tsv` | Deletion 长度分布（按事件长度分组）：Length_bp, Events, Reads, Reads_Pct |
| `insertion_length.tsv` | Insertion 长度分布：Length_bp, Events, Reads, Reads_Pct |
| `event_level_details.tsv` | 事件级统计：每个突变事件一行，记录 Type, Start_Pos, End_Pos, Length, Affected_Targets, N_Targets, Target_Range, Sequences, Reads |

---

## Parameter tuning / 参数调整

| 问题 | 调整方法 |
|---------|-----------|
| 假阳性 indel 过多 | 增大 `--gap-open`（如 -3.0 → -4.0） |
| 遗漏真实 indel | 降低 `--min-scale`（如 1.0 → 0.8） |
| 保守区假阳性点突变 | 增大 `--max-scale`（如 6.0 → 10.0） |
| Indel 被短 match 碎片化 | 开启 `--gap-exit-strength -1.0` |
| 重复序列匹配到错误副本 | 开启 `--homology-penalty -1.0` |
| 纯点突变 allele 假阳性过多 | 增大 `--min-reads-sub`（如 5 → 10） |

### 常用预设

**标准扩增子（单靶点）：**
```
--match-score 2.0 --mismatch-penalty -3.0 --gap-open -2.0 --gap-extend -0.1
```

**多靶点谱系（默认）：**
```
--lineage --min-scale 1.0 --cutsite-edge-scale 2.0 --max-scale 6.0 --sub-window 3
```

**多靶点谱系（严格，更低假阳性）：**
```
--lineage --max-scale 10.0 --sub-window 2 --mismatch-density-threshold 0.40 --min-reads-sub 10
```

---

## Performance notes / 性能说明

- 时间复杂度：标准模式和谱系模式均为 $O(m \times n)$。
- 并行：默认单线程。用 `--threads N` 开启多进程并行。
- 线程数建议不超过 CPU 物理核心数。
- 向量化：Iy 和 M 按行 NumPy 向量化；Ix 因同行数据依赖保持顺序循环。
