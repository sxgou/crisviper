# CRISPR谱系追踪扩增子分析软件 - 详细分析流程说明

## 项目概述

CRISPR谱系追踪扩增子分析软件是一个专门用于分析CRISPR-Cas9 lineage tracing实验中合成靶标阵列的amplicon测序数据的工具。该软件能够准确识别编辑事件，包括大片段缺失、插入和点突变，同时过滤掉测序artifact。

### 核心特性
- **双端reads合并**：支持使用fastp进行双端测序数据合并，可配置最小重叠长度
- **Primer验证**：仅分析包含正确Primer5和Primer3（可配置错配容差）的reads
- **Anchor-guided比对**：使用保守区域作为anchor进行准确比对，即使存在大片段缺失
- **片段分类**：正确分类缺失内的"fragment"序列为MMEJ或NHEJ修复产物
- **位置过滤的SNV检测**：仅在cutsite±3bp窗口内的点突变被视为有效编辑事件
- **靶标状态确定**：为每个靶标生成突变状态矩阵（WT、DEL、INS、SNV、COMPLEX）
- **富集分析**：通过比较断点分布与预期cutsite验证切割特异性

## 整体分析流程

### 1. 数据预处理
- 对于双端测序数据，使用fastp进行reads合并
- 质量过滤和adapter trimming
- 生成合并统计报告

### 2. Primer验证
- 检测read的5'端是否包含Primer5
- 检测read的3'端是否包含Primer3
- 基于错配率阈值判断primers是否有效
- 无效reads被标记为无效，不进行后续分析

### 3. Anchor匹配
- 从配置文件中提取anchor序列（primers、prefix、conserved regions、PAM、postfix）
- 使用k-mer索引快速匹配anchor
- 执行半全局比对确定匹配位置和一致性

### 4. 缺失推断
- 基于matched anchors之间的间隙推断缺失区域
- 计算缺失长度和位置
- 识别与缺失区域对应的read片段

### 5. 突变分类
- 纯缺失：无插入序列
- MMEJ缺失：检测微同源区域
- NHEJ缺失+插入：有插入序列但无微同源
- SNV检测：比对read与参考序列
- SNV过滤：仅保留cutsite窗口内的点突变

### 6. 靶标状态确定
- 为每个靶标分配突变状态
- 基于覆盖度判断靶标是否被删除
- 识别cutsite区域内的断点
- 处理复杂事件（多种突变类型组合）

### 7. Barcode生成
- 将每个靶标的状态编码为字符串
- 生成lineage tracing barcode

## 双端reads合并细节

### 使用的工具：fastp
fastp是一个快速的all-in-one FASTQ预处理器，专门用于处理高通量测序数据。

### 合并参数
```python
# 核心合并函数参数
def merge_paired_end_reads(
    fastq_r1: Path,      # R1 FASTQ文件
    fastq_r2: Path,      # R2 FASTQ文件  
    output_dir: Path,    # 输出目录
    sample_name: str = "sample",
    min_overlap: int = 10,           # 最小重叠长度（bp）
    max_mismatch_rate: float = 0.2,  # 最大错配率（0-1）
    min_quality: int = 20,           # 最低碱基质量
    threads: int = 1,                # 线程数
    verbose: bool = False
)
```

### 合并算法步骤
1. **质量过滤**：移除低质量碱基和reads（基于`--qualified_quality_phred`）
2. **Adapter trimming**：自动检测并移除adapter序列
3. **重叠检测**：在R1和R2之间查找重叠区域
4. **合并**：将重叠的reads合并为单个一致性序列
5. **错误校正**：校正重叠区域的测序错误

### 关键参数说明
| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min_overlap` | 10 | **最小重叠长度（bp）**，低于此值的read pairs不会被合并 |
| `--overlap_diff_limit` | 2 | 重叠区域允许的最大差异数（基于`min_overlap * max_mismatch_rate`计算） |
| `--correction` | 启用 | 启用重叠区域的碱基校正 |
| `--qualified_quality_phred` | 20 | 合格碱基的最低Phred质量分数 |

### 输出文件
- `{sample}_merged.fastq.gz`：合并后的reads
- `{sample}_fastp.json`：JSON格式合并统计报告
- `{sample}_fastp.html`：HTML格式可视化报告

## Primer验证原理

### 数据结构
```python
@dataclass
class PrimerMatch:
    primer_name: str    # "Primer5" 或 "Primer3"
    read_start: int    # read中匹配起始位置
    read_end: int      # read中匹配结束位置
    identity: float    # 匹配一致性（0-1）
    is_valid: bool     # primer是否有效
```

### 检测算法
#### 1. Primer5检测（5'端）
```python
def detect_primer_at_end(
    read_seq: str, 
    primer_seq: str, 
    end: str = 'start',  # 'start'表示5'端
    max_mismatch_rate: float = 0.15,
    min_overlap_ratio: float = 0.8
)
```
- 搜索窗口：`read_seq[:primer_len + 10]`（允许最多10个额外碱基）
- 滑动窗口比对：在搜索窗口内滑动primer长度窗口
- 使用edlib进行快速比对
- 计算编辑距离和一致性

#### 2. Primer3检测（3'端）
- 搜索窗口：`read_seq[-(primer_len + 10):]`
- 搜索偏移量：`max(0, len(read_seq) - len(search_window))`
- 其他步骤与Primer5检测相同

### 验证阈值
- **最大错配率**：`primer_max_mismatch_rate`（默认0.15）
- **最小重叠比例**：`min_overlap_ratio`（默认0.8）
- **最小一致性**：`(primer_len - mismatches) / primer_len`

### 验证结果
```python
def validate_primers(
    read_seq: str,
    primer5_seq: str,
    primer3_seq: str,
    max_mismatch_rate: float = 0.15
) -> Tuple[bool, Optional[Dict], bool, Optional[Dict]]:
```
返回：`(has_primer5, primer5_match_info, has_primer3, primer3_match_info)`

## Anchor-guided比对算法

### 数据结构
```python
@dataclass
class Anchor:
    name: str           # anchor名称
    ref_start: int      # 参考序列起始位置（0-based inclusive）
    ref_end: int        # 参考序列结束位置（0-based exclusive）
    seq: str            # anchor序列
    type: str           # anchor类型：'primer', 'prefix', 'conserved', 'pam', 'postfix'

@dataclass  
class AnchorMatch:
    anchor: Anchor      # 匹配的anchor
    read_start: int     # read中匹配起始位置
    read_end: int       # read中匹配结束位置
    identity: float     # 匹配一致性（0-1）
    cigar: str          # CIGAR字符串
```

### Anchor类型提取
从配置文件features中提取以下anchor类型：
1. **Primers**：Primer5和Primer3
2. **Prefix/Postfix**：靶标阵列前后的保守序列
3. **Conserved regions**：每个靶标的保守区域（13bp）
4. **PAM sequences**：靶标间的PAM_Linker中的PAM区域（3bp）

### AnchorIndex索引构建
```python
class AnchorIndex:
    def __init__(self, k: int = 9):
        self.k = k  # k-mer大小
        self.index: Dict[str, List[Tuple[int, int]]] = {}  # kmer -> [(anchor_id, position)]
        self.anchors: List[Anchor] = []
    
    def build(self, anchors: List[Anchor]):
        # 从每个anchor序列提取所有k-mers
        for anchor_id, anchor in enumerate(anchors):
            anchor_seq = anchor.seq
            for pos in range(len(anchor_seq) - self.k + 1):
                kmer = anchor_seq[pos:pos+self.k]
                if kmer not in self.index:
                    self.index[kmer] = []
                self.index[kmer].append((anchor_id, pos))
```

### Anchor匹配流程
#### 1. k-mer提取和候选anchor识别
```python
candidate_anchors: Set[int] = set()
for i in range(len(read_seq) - k + 1):
    kmer = read_seq[i:i+k]
    matches = anchor_index.query(kmer)  # 查询k-mer索引
    for anchor_id, _ in matches:
        candidate_anchors.add(anchor_id)
```

#### 2. 半全局比对
```python
def semi_global_align(query_seq: str, target_seq: str) -> Dict[str, Any]:
    # 使用edlib的NW（Needleman-Wunsch）模式
    result = edlib.align(query_seq, target_seq, mode="NW", task="path")
    
    # 计算一致性
    edit_distance = result["editDistance"]
    aligned_length = max(len(query_seq), len(target_seq))
    identity = 1.0 - (edit_distance / aligned_length) if aligned_length > 0 else 0.0
    
    return {
        'query_start': query_start,
        'query_end': query_end,
        'identity': identity,
        'edit_distance': edit_distance,
        'cigar': result.get("cigar", "")
    }
```

#### 3. 匹配过滤条件
- **最小一致性**：`min_anchor_identity`（默认0.8）
- **最小对齐长度**：至少anchor长度的70%
- **重叠解决**：保留一致性最高的非重叠匹配

### 缺失推断
```python
def infer_deletions(anchor_matches: List[AnchorMatch], ref_length: int) -> List[Tuple[int, int]]:
    deletions = []
    prev_ref_end = 0
    
    # 按参考位置排序
    sorted_matches = sorted(anchor_matches, key=lambda m: m.anchor.ref_start)
    
    for match in sorted_matches:
        anchor = match.anchor
        if anchor.ref_start > prev_ref_end:
            # anchor之间的间隙即为缺失区域
            deletions.append((prev_ref_end, anchor.ref_start))
        prev_ref_end = max(prev_ref_end, anchor.ref_end)
    
    # 检查末端缺失
    if prev_ref_end < ref_length:
        deletions.append((prev_ref_end, ref_length))
    
    return deletions
```

## 突变检测和分类

### 1. 缺失事件分类

#### 纯缺失检测
```python
if len(read_fragment) == 0:
    return {
        'type': 'pure_deletion',
        'del_length': del_end - del_start,
        'inserted_seq': '',
        'is_mmej': False,
        'mh_left': '',
        'mh_right': ''
    }
```

#### 微同源（MMEJ）检测
```python
def detect_microhomology(fragment: str, flank_left: str, flank_right: str, 
                         min_homology: int = 3, max_mismatch: int = 1):
    # 尝试所有可能的分割点
    for split in range(min_homology, len(fragment) - min_homology + 1):
        left_part = fragment[:split]
        right_part = fragment[split:]
        
        # 检查left_part是否匹配left_flank末端
        if len(flank_left) >= len(left_part):
            left_suffix = flank_left[-len(left_part):]
            mismatches = sum(1 for a, b in zip(left_suffix, left_part) if a != b)
            left_match = mismatches <= max_mismatch
        
        # 检查right_part是否匹配right_flank起始
        if len(flank_right) >= len(right_part):
            right_prefix = flank_right[:len(right_part)]
            mismatches = sum(1 for a, b in zip(right_prefix, right_part) if a != b)
            right_match = mismatches <= max_mismatch
        
        if left_match and right_match:
            return left_part, right_part  # 发现微同源
    
    return None  # 未发现微同源
```

#### 缺失分类结果
- **纯缺失**：无插入序列
- **MMEJ缺失**：检测到微同源，`is_mmej=True`
- **NHEJ缺失+插入**：有插入序列但无微同源，`is_mmej=False`

### 2. 点突变（SNV）检测和过滤

#### SNV检测
```python
def process_point_mutations(read_seq: str, ref_seq: str, anchor_matches: List[Dict], 
                           structure: AmpliconStructure, window_size: int = 3):
    valid_snvs = []
    
    for match in anchor_matches:
        read_start = match.get('read_start', 0)
        ref_start = match.get('ref_start', 0)
        length = match.get('length', 0)
        
        for i in range(length):
            read_pos = read_start + i
            ref_pos = ref_start + i
            
            if read_pos >= len(read_seq) or ref_pos >= len(ref_seq):
                continue
                
            read_base = read_seq[read_pos]
            ref_base = ref_seq[ref_pos]
            
            if read_base != ref_base:
                # 检查是否在有效SNV窗口内
                is_valid, target_name = is_in_valid_snv_window(ref_pos, structure, window_size)
                
                if is_valid:
                    snv = SNVEvent(
                        ref_pos=ref_pos,
                        ref_base=ref_base,
                        alt_base=read_base,
                        target_name=target_name
                    )
                    valid_snvs.append(snv)
    
    return valid_snvs
```

#### 有效SNV窗口检查
```python
def is_in_valid_snv_window(ref_pos: int, structure: AmpliconStructure, window_size: int = 3):
    for target in structure.targets:
        # 计算cutsite周围的窗口
        window_start = max(0, target.cutsite_start - window_size)
        window_end = min(len(structure.reference), target.cutsite_end + window_size)
        
        if window_start <= ref_pos < window_end:
            return True, target.name
    
    return False, None
```

### 3. 靶标状态确定

#### 状态确定算法
```python
def determine_target_state(target: Target, deletion_events: List[DeletionEvent],
                          insertion_events: List[InsertionEvent], snv_events: List[SNVEvent],
                          coverage_threshold: float = 0.9) -> str:
```

#### 状态优先级（从高到低）：
1. **DELETED**：靶标被缺失完全覆盖（覆盖度≥90%）
2. **DEL:{length}**：缺失断点在cutsite区域（±5bp内）
3. **INS:{sequence}**：插入事件在cutsite区域
4. **SNV:{mutation}**：cutsite区域内的点突变
5. **COMPLEX**：多种突变类型组合
6. **WT**：未检测到突变

#### 覆盖度计算
```python
def calculate_coverage(target_start: int, target_end: int, 
                      event_start: int, event_end: int) -> float:
    overlap_start = max(target_start, event_start)
    overlap_end = min(target_end, event_end)
    
    if overlap_end <= overlap_start:
        return 0.0
    
    overlap_length = overlap_end - overlap_start
    target_length = target_end - target_start
    
    return overlap_length / target_length if target_length > 0 else 0.0
```

## 输出文件格式

### 1. per_read_annotation.tsv
每行对应一个read的分析结果，列包括：

| 列名 | 描述 |
|------|------|
| `read_id` | Read标识符 |
| `valid` | Read是否有效（通过primer验证） |
| `Target1` - `Target10` | 每个靶标的状态 |
| `barcode` | lineage tracing barcode（靶标状态用`|`分隔） |
| `events_json` | 突变事件的JSON格式详情 |

**示例**：
```
read_id    valid   Target1  Target2  ...  Target10  barcode                                events_json
read001    True    DELETED  WT       ...  SNV:A>G   DELETED|WT|...|SNV:A>G                {"deletions": [...], "snvs": [...]}
```

### 2. barcode_frequencies.tsv
统计每个unique barcode的出现频率：

| 列名 | 描述 |
|------|------|
| `barcode` | Unique barcode字符串 |
| `count` | 该barcode的出现次数 |
| `frequency` | 相对频率 |

**示例**：
```
barcode                                count  frequency
WT|WT|WT|...|WT                       15000  0.75
DELETED|DELETED|...|DELETED            3000  0.15
WT|DEL:45|WT|...|WT                   1000  0.05
```

### 3. statistics.json
JSON格式的统计摘要：

```json
{
  "total_reads": 20000,
  "passed_reads": 19500,
  "primer5_detected": 19800,
  "primer3_detected": 19700,
  "both_primers_detected": 19500,
  "editing_efficiency": 0.43,
  "per_target_statistics": {
    "Target1": {
      "total": 19500,
      "wt": 11000,
      "deleted": 5000,
      "del_lengths": [45, 67, ...],
      "snv_count": 300,
      "insertion_count": 50
    },
    ...
  },
  "breakpoint_density": [...],
  "indel_length_distribution": [...],
  "snv_distribution": [...]
}
```

## 参数说明

### 核心分析参数

| 参数 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `--primer_max_mismatch_rate` | float | 0.15 | Primer验证最大错配率（0-1） |
| `--min_anchor_identity` | float | 0.8 | Anchor匹配最小一致性 |
| `--kmer_size` | int | 9 | Anchor索引的k-mer大小 |
| `--snv_window` | int | 3 | Cutsite周围的SNV窗口大小（bp） |
| `--min_deletion_size` | int | 20 | 报告为大片段缺失的最小长度（bp） |
| `--trim_primers` | bool | True | 是否在分析前trim primer序列 |

### 双端合并参数

| 参数 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `--min_overlap` | int | 30 | **最小重叠长度（bp）** |
| `--max_mismatch_rate` | float | 0.2 | 重叠区域最大错配率 |
| `--min_quality` | int | 20 | 最低碱基质量阈值 |
| `--skip_merge` | bool | False | 跳过合并（仅分析R1） |
| `--threads` | int | 1 | fastp合并线程数 |

### 输入/输出参数

| 参数 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `--fastq_r1` | str | 必需 | R1 FASTQ文件路径 |
| `--fastq_r2` | str | 可选 | R2 FASTQ文件路径（双端模式） |
| `--config` | str | 必需 | 扩增子结构配置文件 |
| `--output_dir` | str | ./results | 输出目录 |
| `--sample_name` | str | sample | 样本名称 |
| `--max_reads` | int | None | 最大处理reads数（测试用） |

## 配置文件和数据结构

### 配置文件格式（JSON/YAML）
```json
{
  "reference": "ATCG...完整参考序列...",
  "features": [
    {
      "name": "Primer5",
      "start": 0,
      "end": 23,
      "type": "primer"
    },
    {
      "name": "prefix", 
      "start": 23,
      "end": 28,
      "type": "prefix"
    },
    {
      "name": "Target1",
      "start": 28,
      "end": 48,
      "type": "target",
      "conserved": [28, 41],
      "cutsite": [41, 48]
    },
    {
      "name": "PAM_Linker1",
      "start": 48,
      "end": 55,
      "type": "pam_linker",
      "pam": [48, 51],
      "linker": [51, 55]
    },
    ...
    {
      "name": "Primer3",
      "start": 299,
      "end": 332,
      "type": "primer"
    }
  ]
}
```

### 标准扩增子结构
基于典型CRISPR lineage tracing实验设计：

| 组件 | 长度（bp） | 描述 |
|------|------------|------|
| Primer5 | 23 | 5'端primer |
| Prefix | 5 | 5'端保守序列 |
| Target1-10 | 每个20bp | 10个靶标，每个包含：<br>- Conserved region: 13bp<br>- Cutsite region: 7bp（包含上游3bp PAM） |
| PAM_Linker1-9 | 每个7bp | 靶标间的连接序列：<br>- PAM: 3bp<br>- Linker: 4bp |
| Postfix | 8 | 3'端保守序列 |
| Primer3 | 33 | 3'端primer |

**总长度**：332bp

## 性能优化和注意事项

### 1. 内存使用优化
- **流式处理**：逐个read处理，避免同时加载所有reads到内存
- **k-mer索引**：使用高效哈希表存储anchor k-mers
- **CIGAR压缩**：仅存储必要的比对信息

### 2. 计算效率优化
- **edlib集成**：使用C++编写的edlib库进行快速序列比对
- **多线程支持**：fastp合并支持多线程
- **向量化操作**：numpy加速数值计算

### 3. 算法精度保证
- **锚点验证**：多个anchor的协同验证提高比对准确性
- **微同源检测**：动态规划检测MMEJ事件
- **位置过滤**：严格的SNV窗口过滤减少假阳性

### 4. 常见问题解决

#### 缺失检测错误
- **问题**：大片段缺失导致anchor匹配失败
- **解决**：使用多个保守区域作为anchor，即使部分anchor缺失也能比对

#### Primer错配
- **问题**：测序错误导致primer检测失败
- **解决**：可配置错配容差（`--primer_max_mismatch_rate`）

#### 重叠reads合并
- **问题**：插入片段长度变异影响重叠检测
- **解决**：调整`--min_overlap`参数适应不同插入片段分布

## 扩展性和自定义

### 添加新anchor类型
在`_extract_anchors()`方法中添加新的anchor类型检测逻辑：

```python
if feat_type == 'new_feature_type':
    seq = self.structure.reference[start:end]
    anchor = Anchor(
        name=name,
        ref_start=start,
        ref_end=end,
        seq=seq,
        type='new_type'
    )
    anchors.append(anchor)
```

### 自定义突变分类
在`mutation_classification.py`中添加新的分类算法：

```python
def classify_custom_event(read_fragment: str, ref_seq: str, del_start: int, del_end: int):
    # 实现自定义分类逻辑
    return {
        'type': 'custom_event',
        'custom_field': 'value'
    }
```

### 修改输出格式
修改`engine.py`中的`analyze_fastq()`方法以生成自定义输出格式。

## 参考文献和算法来源

### 核心算法参考
1. **fastp**：Shifu Chen, et al. (2018). "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics.
2. **edlib**：M. Šošić, M. Šikić (2017). "Edlib: a C/C++ library for fast, exact sequence alignment using edit distance."
3. **MMEJ检测**：基于微同源介导的末端连接（Microhomology-Mediated End Joining）机制

### 生物学背景
- **CRISPR-Cas9 lineage tracing**：通过合成靶标阵列追踪细胞谱系
- **DNA修复机制**：NHEJ、MMEJ等DNA双链断裂修复途径
- **扩增子测序**：目标区域扩增后的高通量测序

### 质量控制指标
- **编辑效率**：`valid_editing_reads / total_valid_reads`
- **断点富集**：cutsite区域断点密度与背景比较
- **突变分布**：缺失长度、SNV位置的分布分析

---

## 附录：命令行使用示例

### 完整分析流程
```bash
# 双端数据分析（合并+分析）
lineage-tracer analyze \
  --fastq_r1 sample_R1.fastq.gz \
  --fastq_r2 sample_R2.fastq.gz \
  --config amplicon_structure.json \
  --output_dir ./results \
  --sample_name Embryo_Day5 \
  --min_overlap 30 \
  --max_mismatch_rate 0.2 \
  --threads 8

# 单端数据分析
lineage-tracer analyze \
  --fastq_r1 merged.fastq.gz \
  --config amplicon_structure.json \
  --output_dir ./results \
  --sample_name Merged_Sample \
  --primer_max_mismatch_rate 0.15

# 配置向导
lineage-tracer config-wizard --output my_config.json
```

### 结果验证
```bash
# 查看统计摘要
cat results/statistics.json | python -m json.tool

# 查看前10个reads的注释
head -11 results/per_read_annotation.tsv | column -t -s $'\\t'

# 查看barcode频率
sort -k3 -nr results/barcode_frequencies.tsv | head -20
```

---
*文档版本：1.0.0*
*最后更新：2024年4月21日*
*项目版本：lineage_tracer_amplicon v1.0.0*
