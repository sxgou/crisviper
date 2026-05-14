# Lineage Tracer — 谱系示踪序列分析工具

[![CI](https://github.com/NousResearch/lineage-tracer/actions/workflows/ci.yml/badge.svg)](https://github.com/NousResearch/lineage-tracer/actions/workflows/ci.yml)

基于仿射gap惩罚算法的序列比对命令行工具，支持**结构感知的谱系示踪比对模式**、多进程并行加速和FASTQ格式转换。专为多靶点基因编辑谱系示踪实验设计。

## 功能特性

- **标准序列比对**：基于Gotoh算法（Needleman-Wunsch的仿射gap扩展）
- **谱系示踪比对模式**：结构感知的靶向分析，位置依赖gap惩罚，优先在cutsite区域开启gap
- **假阳性突变过滤**：自动矫正cutsite区域外的点突变，保留紧邻gap的突变
- **高密度mismatch检测**：>34% mismatch密度的区域自动转换为indel
- **跨靶点重复序列矫正**：自动检测并修正因重复序列导致的错误比对（11组重复序列，4-15bp）
- **小片段跨靶点矫正**：矫正TAGTAT和单碱基A的跨靶点错误匹配
- **孤立匹配清除**：合并被偶然匹配碱基打断的连续deletion片段
- **多进程并行加速**：自动利用所有CPU核心，大幅提升批量比对速度
- **FASTQ格式转换**：支持将FASTQ文件转换为TSV或FASTA格式
- **灵活的比对参数**：可调整匹配得分、错配惩罚、gap惩罚等参数
- **多种输出格式**：支持JSON和TSV输出格式
- **分析报告生成**：支持JSON和HTML格式的突变分析报告
- **实时进度显示**：并行处理时实时显示完成进度

## 快速开始

### 安装依赖
```bash
pip install -r requirements.txt
```

### 基本使用示例

1. **将FASTQ转换为TSV格式**：
```bash
python carlin_tool.py convert fastq-to-tsv \
  --fastq reads.fastq.gz \
  --output reads.tsv \
  --sample-name my_sample
```

2. **标准并行批量序列比对**（自动使用所有CPU核心）：
```bash
python carlin_tool.py align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json
```

3. **谱系示踪比对模式**（结构感知，自动推断cutsite位置）：
```bash
python carlin_tool.py align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --lineage
```

4. **指定并行进程数**：
```bash
python carlin_tool.py align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json \
  --threads 8
```

## 文件结构

```
.
├── affine_gap_alignment.py    # 核心比对算法（含标准算法和谱系示踪算法，跨靶点矫正）
├── carlin_tool.py             # 命令行工具主程序
├── run_corrected.py           # 批量比对运行器（谱系模式，含HTML报告生成）
├── requirements.txt           # Python依赖
├── README.md                  # 说明文档
├── example_data/              # 示例数据
│   ├── test.fastq.gz         # FASTQ测序数据示例
│   └── test_queries.tsv      # 查询序列示例
├── results/                   # 输出结果目录
│   ├── correct.json          # 矫正后的比对结果
│   └── report.html           # HTML分析报告
└── docs/                      # 详细文档目录
    ├── ALGORITHM.md          # 算法原理文档（含跨靶点矫正算法）
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

## 并行加速

- 使用Python多进程（`ProcessPoolExecutor`）实现真正的并行计算
- 每个CPU核心独立处理不同的序列，互不干扰
- 默认自动检测并使用所有可用CPU核心
- 可通过 `--threads` 或 `-t` 参数手动指定进程数
- 支持对数千条序列的快速批量处理

## 应用场景

- **谱系示踪分析**：多靶点基因编辑的突变组合检测
- **CARLIN基因编辑分析**：检测大片段缺失和点突变
- **扩增子测序分析**：比对测序reads到参考序列
- **突变检测**：识别序列中的插入、缺失和替换
- **序列质量控制**：评估测序reads与参考序列的相似度

## 性能特点

- **时间复杂度**：O(m×n)，其中m,n为序列长度
- **并行加速比**：近似线性（受序列长度和数据量影响）
- **典型性能**：332bp×133bp（CARLIN）单序列约45ms

## 许可证

MIT License

## 技术支持

如有问题或建议，请：
1. 确保已安装所有依赖 (`pip install -r requirements.txt`)
2. 查看详细文档：`docs/ALGORITHM.md` 和 `docs/USER_GUIDE.md`
3. 使用帮助命令：`python carlin_tool.py --help`
4. 谱系示踪模式：`python carlin_tool.py align --help`
