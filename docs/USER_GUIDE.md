# User guide / 使用手册

## Installation / 安装

```bash
pip install crisviper
```

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

`--lineage` 开启结构感知的 gap 惩罚、自动检测 cutsite 位置、后处理矫正。HTML 报告包含突变类型分布、长度分布和摘要统计。

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
  --gap-exit-bonus -1.0 \
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
| `--fastq` | — | 输入 FASTQ（支持 .gz） |
| `--output` | — | 输出 TSV 路径 |
| `--sample-name` | `sample` | 样本名称 |
| `--min-reads` | 1 | 最小 read 数过滤 |

### `crisviper convert fastq-to-fasta`

同上，但不支持 `--min-reads`。

### `crisviper align`

**基本选项：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--reference` | — | 参考序列 FASTA |
| `--queries` | — | 查询序列文件（TSV 或 FASTA） |
| `--output` | — | 输出路径 |
| `--format` | `json` | 输出格式：`json`、`tsv`、`all` |
| `--threads` / `-t` | 1 | 并行进程数（默认 1=单线程） |

**打分参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--match-score` | 2.0 | 匹配得分 |
| `--mismatch-penalty` | -3.0 | 错配惩罚 |
| `--gap-open` | -2.0 | Gap 开启惩罚 |
| `--gap-extend` | -0.1 | Gap 延伸惩罚 |
| `--global` | off | 全局比对（默认半全局） |

**谱系模式参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--lineage` | off | 开启谱系示踪比对 |
| `--cutsite-scale` | 1.0 | Cutsite gap 惩罚倍率（越低越易开 gap） |
| `--flank-scale` | 2.0 | 侧翼区倍率（±flank-width） |
| `--far-scale` | 6.0 | 保守区倍率 |
| `--flank-width` | 3 | 侧翼宽度 (bp) |
| `--mutation-window` | 3 | Cutsite 保留点突变的窗口半径 (bp) |
| `--density-threshold` | 0.34 | Mismatch 密度阈值，超过则转为 indel |
| `--cutsites` | auto | Cutsite JSON 文件路径（省略则自动检测） |

**DP 原生特征（谱系模式推荐）：**

| Option | Default | Recommended | Description / 说明 |
|--------|---------|-------------|-------------|
| `--gap-exit-bonus` | 0.0 | -1.0 | Gap→match 转换的惩罚；合并碎片化 indel |
| `--short-match-window` | 0 | 3 | 短匹配区域阈值 (bp)，0=关闭 |
| `--short-match-discount` | 1.0 | 0.5 | 短匹配区域得分折扣 (1.0=不打折) |
| `--dense-mismatch-window` | 6 | 6 | 密集错配检测窗口 (bp) |
| `--dense-mismatch-penalty` | 0.0 | -2.0 | 密集错配区域额外惩罚，0=关闭 |
| `--homology-window` | 8 | 8 | 同源区域检测窗口 (bp) |
| `--homology-penalty` | 0.0 | -1.0 | 同源区域惩罚，0=关闭 |
| `--isolated-base-penalty` | 0.0 | -2.0 | 孤立碱基匹配惩罚，0=关闭 |

**矫正管线控制：**

| Option | Effect / 效果 |
|--------|--------|
| `--repeat-correction-mode` | `auto`（默认）、`hardcoded` 或 `off` |
| `--disable-target-misalignment` | 关闭 TAGTAT / 单碱基 A 的跨靶点矫正 |
| `--disable-isolated-match-removal` | 保留切割 deletion 的孤立匹配 |
| `--disable-dense-mismatch-correction` | 关闭后处理密集错配→indel 转换 |
| `--disable-point-mutation-filtering` | 保留所有点突变，不按位置过滤 |

**引物参数：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--primer5-len` | 23 | 5' 端引物长度 (bp) |
| `--primer3-len` | 33 | 3' 端引物长度 (bp) |
| `--primer5-threshold` | 19 | 5' 引物匹配碱基数阈值 |
| `--primer3-threshold` | 29 | 3' 引物匹配碱基数阈值 |

**报告选项：**

| Option | Default | Description / 说明 |
|--------|---------|-------------|
| `--report` | — | 生成报告：`json` 或 `html` |
| `--report-output` | auto | 报告路径（默认从 `--output` 推断） |
| `--min-reads` | 1 | allele 过滤的最小 read 数 |
| `--allele-top-n` | 50 | 报告中展示的 top N alleles |
| `--allele-window-start` | 0 | Allele 热图显示起始位置 |
| `--allele-window-end` | auto | Allele 热图显示结束位置（含） |

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
    "n_mutations_corrected": 1,
    "dense_regions_converted": false
  }
}
```

### TSV

简化的表格：readName, cellBC, UMI, readCount, score, matches, mismatches, gaps_in_query, similarity, aligned_ref, aligned_query, error。

### HTML report / 分析报告

`--report html` 生成的 HTML 报告包含：摘要统计、突变类型条形图、长度分布、top allele 表和每条突变序列的突变标签。所有 JavaScript 内联，无外部依赖，离线可用。

---

## Parameter tuning / 参数调整

| 问题 | 调整方法 |
|---------|-----------|
| 假阳性 indel 过多 | 增大 `--gap-open`（如 -3.0 → -4.0） |
| 遗漏真实 indel | 降低 `--cutsite-scale`（如 1.0 → 0.8） |
| 保守区假阳性点突变 | 增大 `--far-scale`（如 6.0 → 10.0） |
| Indel 被短 match 碎片化 | 开启 `--gap-exit-bonus -1.0` |
| 重复序列匹配到错误副本 | 开启 `--homology-penalty -1.0` |

### 常用预设

**标准扩增子（单靶点）：**
```
--match-score 2.0 --mismatch-penalty -3.0 --gap-open -2.0 --gap-extend -0.1
```

**多靶点谱系（默认）：**
```
--lineage --cutsite-scale 1.0 --flank-scale 2.0 --far-scale 6.0 --mutation-window 3
```

**多靶点谱系（严格，更低假阳性）：**
```
--lineage --far-scale 10.0 --mutation-window 2 --density-threshold 0.40
```

---

## Performance notes / 性能说明

- 时间复杂度：标准模式和谱系模式均为 $O(m \times n)$。
- 并行：默认单线程。用 `--threads N` 开启多进程并行。
- 线程数建议不超过 CPU 物理核心数。
- 向量化：Iy 和 M 按行 NumPy 向量化；Ix 因同行数据依赖保持顺序循环。
