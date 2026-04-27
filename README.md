# CARLIN序列分析工具

基于仿射gap惩罚算法的序列比对命令行工具，支持多进程并行加速和FASTQ格式转换。

## 功能特性

- **仿射gap惩罚算法**：基于Gotoh算法（Needleman-Wunsch的仿射gap扩展）
- **多进程并行加速**：自动利用所有CPU核心，大幅提升批量比对速度
- **FASTQ格式转换**：支持将FASTQ文件转换为TSV或FASTA格式
- **灵活的比对参数**：可调整匹配得分、错配惩罚、gap惩罚等参数
- **多种输出格式**：支持JSON和TSV输出格式
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

2. **并行批量序列比对**（自动使用所有CPU核心）：
```bash
python carlin_tool.py align \
  --reference reference.fasta \
  --queries queries.tsv \
  --output alignments.json
```

3. **指定并行进程数**：
```bash
python carlin_tool.py align \
  --reference reference.fasta \
  --queries queries.fasta \
  --output alignments.json \
  --threads 8
```

## 文件结构

```
.
├── affine_gap_alignment.py    # 核心比对算法
├── carlin_tool.py             # 命令行工具主程序
├── requirements.txt           # Python依赖
├── README.md                  # 说明文档
├── example_data/              # 示例数据
│   ├── reference.fa          # 参考序列示例
│   └── test_queries.tsv      # 查询序列示例
└── docs/                      # 详细文档目录
    ├── ALGORITHM.md          # 算法原理文档
    └── USER_GUIDE.md         # 使用手册
```

## 算法特点

1. **最小化错配**：错配惩罚高于gap惩罚，优先选择插入/缺失而非错配
2. **减少短gap**：gap延伸惩罚极低，鼓励连续长gap而非多个短gap
3. **鼓励插入**：gap开启惩罚相对较低，使插入/缺失比错配更有利
4. **处理连续indel**：仿射gap模型自然处理连续插入/缺失事件
5. **支持半全局比对**：允许序列两端自由gap（不惩罚）

## 并行加速

- 使用Python多进程（`ProcessPoolExecutor`）实现真正的并行计算
- 每个CPU核心独立处理不同的序列，互不干扰
- 默认自动检测并使用所有可用CPU核心
- 可通过 `--threads` 或 `-t` 参数手动指定进程数
- 支持对数千条序列的快速批量处理

## 应用场景

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
