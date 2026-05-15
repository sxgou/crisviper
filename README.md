# Lineage Tracer — 谱系示踪序列分析工具

[![CI](https://github.com/NousResearch/lineage-tracer/actions/workflows/ci.yml/badge.svg)](https://github.com/NousResearch/lineage-tracer/actions/workflows/ci.yml)

基于仿射gap惩罚算法的序列比对命令行工具，支持**结构感知的谱系示踪比对模式**、**NumPy向量化DP加速**、多进程并行加速和FASTQ格式转换。专为多靶点基因编辑谱系示踪实验设计。

## 功能特性

### 核心算法
- **标准序列比对**：基于Gotoh算法（Needleman-Wunsch的仿射gap扩展）
- **谱系示踪比对模式**：结构感知的靶向分析，位置依赖gap惩罚，优先在cutsite区域开启gap
- **向量化DP加速**：NumPy按行向量化Iy/M递推，性能提升4.3x

### DP原生特征（谱系模式）
- **Gap Exit Bonus**：gap→match转换时额外奖励，抑制碎片化indel
- **Short Match Discount**：短匹配区域降低match奖励，吸收小片段到gap中
- **Dense Mismatch Penalty**：连续密集错配区域附加惩罚，驱动indel形成
- **Homology Penalty**：同源区域降低match奖励，抑制重复序列错误匹配
- **Isolated Base Penalty**：孤立匹配碱基附加惩罚，避免碱基分割deletion块
- **Cumsum密度计算**：基于对角线的cumsum代替O(m×n×window)三重循环

### 后处理矫正管线
- **假阳性突变过滤**：自动矫正cutsite区域外的点突变，保留紧邻gap的突变
- **高密度mismatch检测**：>34% mismatch密度的区域自动转换为indel
- **跨靶点重复序列矫正**：自动检测并修正因重复序列导致的错误比对（11组重复序列，4-15bp）
- **小片段跨靶点矫正**：矫正TAGTAT和单碱基A的跨靶点错误匹配
- **孤立匹配清除**：合并被偶然匹配碱基打断的连续deletion片段

### 工程特性
- **多进程并行加速**：自动利用CPU核心并行处理（默认上限12线程防止系统过载）
- **ProcessPoolExecutor稳定性**：OMP_NUM_THREADS=1 + BrokenProcessPool自动单线程回退
- **FASTQ格式转换**：支持将FASTQ文件转换为TSV或FASTA格式
- **灵活的比对参数**：可调整匹配得分、错配惩罚、gap惩罚等参数
- **多种输出格式**：支持JSON和TSV输出格式
- **分析报告生成**：支持JSON和HTML格式的突变分析报告（含可视化图表、突变标签）
- **实时进度显示**：并行处理时实时显示完成进度

## 快速开始

### 安装依赖
```bash
pip install -r requirements.txt
```

### 基本使用示例

1. **将FASTQ转换为TSV格式**：
```bash
lineage-tracer convert fastq-to-tsv \
  --fastq reads.fastq.gz \
  --output reads.tsv \
  --sample-name my_sample
```

2. **标准并行批量序列比对**（自动使用所有CPU核心）：
```bash
lineage-tracer align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json
```

3. **谱系示踪比对模式**（结构感知，自动推断cutsite位置）：
```bash
lineage-tracer align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --lineage
```

4. **谱系模式 + DP原生特征 + 报告**：
```bash
lineage-tracer align \
  --reference example_data/reference.fa \
  --queries example_data/test_queries.tsv \
  --output results/prefix \
  --format all \
  --lineage \
  --gap-exit-bonus -1.0 \
  --short-match-window 3 --short-match-discount 0.5 \
  --dense-mismatch-penalty -2.0 \
  --homology-penalty -1.0 \
  --isolated-base-penalty -2.0 \
  --report html
```

5. **指定并行进程数**：
```bash
lineage-tracer align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --threads 8
```

## 文件结构

```
.
├── ltlib/                     # 核心库（模块化包结构）
│   ├── cli.py                 # 命令行工具主程序（lineage-tracer命令入口）
│   ├── alignment.py           # DP比对算法（标准+位置感知+向量化）
│   ├── lineage.py             # 谱系示踪：gap profile, cutsite检测
│   ├── pipeline.py            # 管线编排（8步流程+ProcessPoolExecutor）
│   ├── corrections.py         # 后处理矫正（重复序列/跨靶点/孤立匹配）
│   ├── mutation.py            # 突变识别与分类
│   ├── models.py              # 类型安全的数据模型
│   ├── config.py              # 数据结构定义
│   ├── io.py                  # FASTQ/TSV/FASTA I/O
│   ├── reporting.py           # 报告生成（JSON+HTML+TSV）
│   ├── plotting.py            # 可视化图表
│   ├── metrics.py             # 多样性指标
│   ├── denoiser.py            # 降噪
│   ├── caller.py              # 等位基因调用
│   ├── threshold.py           # 阈值计算
│   └── logging_config.py      # 日志框架
├── requirements.txt           # Python依赖
├── README.md                  # 说明文档
├── pyproject.toml             # 包构建配置
├── Dockerfile                 # Docker镜像构建
├── config.yaml                # 配置文件模板
├── example_data/              # 示例数据
│   ├── reference.fa          # CARLIN参考序列
│   ├── test_queries.tsv      # 查询序列示例（23,430条）
│   └── test.fastq.gz         # FASTQ测序数据示例
├── results/                   # 测试结果
│   ├── run_comparison.py     # 对比测试脚本
│   ├── comparison.json       # 矫正ON/OFF对比摘要
│   ├── corrected/            # 矫正开启结果
│   │   ├── result.json/tsv  # 比对结果
│   │   └── report.html      # HTML分析报告
│   └── uncorrected/          # 矫正关闭结果
│       ├── result.json/tsv
│       └── report.html
├── tests/                     # 测试套件（199+项）
│   ├── test_gap_exit_bonus.py
│   └── ...
├── workflow/                  # CI/CD工作流
└── docs/                      # 详细文档
    ├── ALGORITHM.md          # 算法原理文档
    └── USER_GUIDE.md         # 使用手册
```

## 算法特点

### 标准比对模式
1. **最小化错配**：错配惩罚高于gap惩罚，优先选择插入/缺失而非错配
2. **减少短gap**：gap延伸惩罚极低，鼓励连续长gap而非多个短gap
3. **处理连续indel**：仿射gap模型自然处理连续插入/缺失事件
4. **支持半全局比对**：允许序列两端自由gap（不惩罚）

### 谱系示踪比对模式（`--lineage`）
1. **结构感知gap惩罚**：cutsite区域gap惩罚降低、保守区域提升，gap自动优先在cutsite开启
2. **突变区域过滤**：cutsite ±3bp 范围外的点突变自动矫正为参考序列
3. **高密度mismatch转换**：连续区域>34%为mismatch的，整段转换为insertion事件
4. **gap邻域保护**：紧邻gap的突变不受矫正（避免假阴性）
5. **自动结构识别**：自动检测10个target的cutsite位置（标准CARLIN扩增子）
6. **跨靶点重复序列矫正**：检测并修正11组因重复序列导致的错误比对（15bp-4bp）
7. **小片段跨靶点矫正**：TAGTAT和单碱基A的跨靶点搬迁
8. **孤立匹配清除**：合并被偶然匹配打断的连续deletion片段

### DP原生特征（v3.0+）
| 特征 | 描述 | 效果 |
|------|------|------|
| Gap Exit Bonus | gap→match转换时附加奖励 | 抑制碎片化indel，倾向合并gap |
| Short Match Discount | 短匹配（<3bp）降低match奖励 | 吸收小片段到gap中 |
| Dense Mismatch Penalty | 密集错配区域附加惩罚 | 驱动DP产生indel而非多错配 |
| Homology Penalty | 同源区域降低match奖励 | 抑制重复序列的过度匹配 |
| Isolated Base Penalty | 孤立匹配碱基附加惩罚 | 避免单碱基分割连续deletion |

## 并行加速

- 使用Python多进程（`ProcessPoolExecutor`）实现真正的并行计算
- 每个CPU核心独立处理不同的序列，互不干扰
- 默认自动检测并使用所有可用CPU核心（上限12线程，防止系统过载）
- 可通过 `--threads` 或 `-t` 参数手动指定进程数
- 支持对数千条序列的快速批量处理

## 性能特点

- **时间复杂度**：O(m×n)，其中m,n为序列长度
- **向量化加速**：NumPy Iy/M按行向量化 + cumsum对角线密度计算
- **典型性能**（332bp CARLIN，全部特性开启）：
  - 单序列：~0.37s（优化前1.6s，4.3x提升）
  - 12线程批处理500条：~19s
  - 12线程全量23,430条：~14min

## 应用场景

- **谱系示踪分析**：多靶点基因编辑的突变组合检测
- **CARLIN基因编辑分析**：检测大片段缺失和点突变
- **扩增子测序分析**：比对测序reads到参考序列
- **突变检测**：识别序列中的插入、缺失和替换
- **序列质量控制**：评估测序reads与参考序列的相似度

## 许可证

MIT License

## 技术支持

如有问题或建议，请：
1. 确保已安装所有依赖 (`pip install -r requirements.txt`)
2. 查看详细文档：`docs/ALGORITHM.md` 和 `docs/USER_GUIDE.md`
3. 使用帮助命令：`lineage-tracer --help`
4. 谱系示踪模式：`lineage-tracer align --help`
