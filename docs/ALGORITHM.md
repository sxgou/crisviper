# 仿射gap惩罚序列比对算法原理

## 概述

本工具实现了两种比对算法：

1. **标准仿射gap比对**：基于Gotoh算法的仿射gap惩罚序列比对，Needleman-Wunsch的仿射gap扩展。
2. **谱系示踪比对**：结构感知的改进算法，支持位置依赖的gap惩罚、高密度mismatch转换和区域感知突变过滤。

专为生物信息学序列比对设计，特别适用于CARLIN等多靶点谱系示踪实验的序列分析。

---

# 第一部分：标准仿射gap比对

## 算法特点

### 1. 仿射gap惩罚模型
- **Gap opening惩罚** ($g_o$): 开启新gap的成本
- **Gap extension惩罚** ($g_e$): 延伸现有gap的成本
- **总gap成本** = $g_o + k \cdot g_e$，其中$k$是gap长度

### 2. 设计原则
1. **最小化错配**：错配惩罚高于gap惩罚，优先选择插入/缺失而非错配
2. **减少短gap**：gap延伸惩罚极低，鼓励连续长gap而非多个短gap
3. **鼓励插入**：gap开启惩罚相对较低，使插入/缺失比错配更有利
4. **处理连续indel**：仿射gap模型自然处理连续插入/缺失事件

## 数学原理

### 动态规划递推

对于序列$X$（长度$m$）和$Y$（长度$n$），定义三个DP矩阵：

#### 1. 匹配/错配状态矩阵 $M(i,j)$
以匹配或错配结束的最佳得分：
$$
M(i,j) = s(x_i, y_j) + \max\begin{cases}
M(i-1, j-1) \\
I_x(i-1, j-1) \\
I_y(i-1, j-1)
\end{cases}
$$

其中 $s(x_i, y_j)$ 是位置得分：
$$
s(x_i, y_j) = 
\begin{cases}
\text{match\_score} & \text{if } x_i = y_j \\
\text{mismatch\_penalty} & \text{otherwise}
\end{cases}
$$

#### 2. X中gap状态矩阵 $I_x(i,j)$
在序列$X$中有gap（$Y$中插入）结束的最佳得分：
$$
I_x(i,j) = \max\begin{cases}
M(i, j-1) + g_o + g_e \\
I_x(i, j-1) + g_e \\
I_y(i, j-1) + g_o + g_e
\end{cases}
$$

#### 3. Y中gap状态矩阵 $I_y(i,j)$
在序列$Y$中有gap（$X$中缺失）结束的最佳得分：
$$
I_y(i,j) = \max\begin{cases}
M(i-1, j) + g_o + g_e \\
I_x(i-1, j) + g_o + g_e \\
I_y(i-1, j) + g_e
\end{cases}
$$

### 初始化条件

#### 半全局比对（默认）
- $M(0,0) = 0$
- $I_x(0,0) = I_y(0,0) = -\infty$
- 第一行：$M(0,j) = 0$, $I_x(0,j) = 0$, $I_y(0,j) = -\infty$
- 第一列：$M(i,0) = 0$, $I_x(i,0) = -\infty$, $I_y(i,0) = 0$

#### 全局比对
- $M(0,0) = 0$
- $I_x(0,0) = I_y(0,0) = g_o$
- 第一行：$M(0,j) = -\infty$, $I_x(0,j) = g_o + (j-1)g_e$, $I_y(0,j) = -\infty$
- 第一列：$M(i,0) = -\infty$, $I_x(i,0) = -\infty$, $I_y(i,0) = g_o + (i-1)g_e$

### 回溯算法

从最终得分最高的位置开始，根据状态转移反向追踪：

1. **匹配状态 (M)**：
   - 添加 $x_i$ 和 $y_j$ 到比对结果
   - 移动到 $(i-1, j-1)$，选择前一个得分最高的状态

2. **X中gap状态 (I_x)**：
   - 添加 '-' 和 $y_j$ 到比对结果
   - 移动到 $(i, j-1)$，选择前一个得分最高的状态

3. **Y中gap状态 (I_y)**：
   - 添加 $x_i$ 和 '-' 到比对结果
   - 移动到 $(i-1, j)$，选择前一个得分最高的状态

---

# 第二部分：谱系示踪比对

## 设计动机

多靶点谱系示踪实验（如CARLIN）使用10个串联的target位点，每个target包含13bp保守区和7bp cutsite。当sgRNA诱导切割后，DNA修复会产生以下突变类型：

- **大片段缺失**：一个或多个target被删除
- **小indel**：cutsite处的插入或缺失（1-20bp）
- **点突变**：主要出现在cutsite附近

生物学先验知识：
1. Gap（indel）主要在cutsite处发生，而非保守/骨架区域
2. 真实点突变也同样集中在cutsite附近
3. 远离cutsite的mismatch很可能是测序错误或PCR错误
4. 连续的高密度mismatch区域实际上是insertion事件

标准Gotoh算法使用**全局固定gap惩罚**，无法利用这些先验知识。

## 位置感知的仿射gap惩罚

### 核心思想

gap开启惩罚不再全局固定，而是**依赖于gap在参考序列上的位置**：

$$
g_o(p) = g_o^{\text{base}} \times s(p), \quad g_e(p) = g_e^{\text{base}} \times s(p)
$$

其中 $s(p)$ 是位置 $p$ 的惩罚倍率：

| 区域 | 倍率 $s(p)$ | 有效gap_open | 含义 |
|------|-------------|-------------|------|
| Cutsite内部 | 1.0 | -2.0 | 鼓励gap检测 |
| 侧翼 ±3bp | 2.0 | -4.0 | 中等抑制 |
| 保守/骨架区 | 2.0（默认） | -4.0 | 一定抑制（用户可通过 `far_scale` 参数调节，增大则增强压制） |

### 修改后的DP递推

在 $I_x$ 和 $I_y$ 状态转移中，使用位置依赖的惩罚：

$$
I_x(i,j) = \max\begin{cases}
M(i, j-1) + g_o(i) + g_e(i) \\
I_x(i, j-1) + g_e(i) \\
I_y(i, j-1) + g_o(i) + g_e(i)
\end{cases}
$$

$$
I_y(i,j) = \max\begin{cases}
M(i-1, j) + g_o(i-1) + g_e(i-1) \\
I_x(i-1, j) + g_o(i-1) + g_e(i-1) \\
I_y(i-1, j) + g_e(i-1)
\end{cases}
$$

其中 $g_o(k)$ 和 $g_e(k)$ 是参考序列第 $k$ 个碱基处的gap惩罚。

### 惩罚倍率的空间分布

```python
# cutsite处: 基准惩罚 1x
for pos in range(cs_start, cs_end + 1):
    gap_open[pos] = base_gap_open * 1.0

# 侧翼 ±flank_width: 2x 惩罚
for offset in range(1, flank_width + 1):
    gap_open[cs_start - offset] = base_gap_open * 2.0
    gap_open[cs_end + offset] = base_gap_open * 2.0

# 其余位置: far_scale 惩罚 (默认 2x)
gap_open[:] = base_gap_open * far_scale  # 初始化

# 然后覆盖cutsite和侧翼的较小值
```

## 高密度mismatch → indel 转换

### 原理

当一段连续比对区域中超过34%的碱基为mismatch时，生物学上更可能是**query序列在该位置有一段insertion**，而非多个独立的点突变。

### 检测算法

使用滑动窗口（默认6bp）扫描alignment：

1. 计算窗口内mismatch密度：$d = \frac{n_{\text{mismatch}}}{n_{\text{bases}}}$
2. 如果 $d > 0.34$ 且 $n_{\text{mismatch}} \geq 2$，标记为密集区
3. 将密集区内所有ref碱基替换为gap（`-`）
4. query碱基保持不变

### 设计决策

- 使用**严格大于**阈值（`> 0.34` 而非 `>= 0.34`），阈值略高于1/3以排除边界情况
- 要求至少2个mismatch，避免孤立点突变被误转换
- 窗口6bp：3/6 = 50% > 34%，2/6 = 33% < 34%，确保足够的mismatch密度

## 区域感知点突变过滤

### 过滤规则

根据谱系示踪实验的生物学特性：

```
                    假阳性矫正区
    ┌──────────────────┐  ┌──────────────────┐
    │  保守区 (矫正)    │  │  保守区 (矫正)    │
    └──────────────────┘  └──────────────────┘
         ┌─┬─┬─┬─┬─┬─┬─┐      ...每个target
         │C│u│t│s│i│t│e│  cutsite 7bp
         └─┴─┴─┴─┴─┴─┴─┘
         ←── ±3bp窗口 ──→  保留点突变
```

**保留条件**（满足任一即保留）：
1. 突变位置在任意cutsite ±3bp范围内
2. 突变位置紧邻一个gap（左侧或右侧相邻位置为gap）

**不满足以上条件** → 将query碱基矫正为ref碱基

### 例外规则：gap邻域保护

紧邻gap的点突变不被矫正，原因：
- Indel边界处的碱基在修复过程中可能同时发生替换
- 测序时indel附近的碱基质量可能较低，但不一定是错误
- 避免将indel-associated SNP误判为假阳性

## 扩增子结构推断

标准CARLIN扩增子结构（332bp）：

```
Primer5(23bp) + prefix(5bp) + 
  [Target(13bp conserved + 7bp cutsite) + PAM_Linker(7bp)] × 9 +
  Target10(13bp + 7bp) + postfix(8bp) + Primer3(33bp)
```

合计：23 + 5 + 10×20 + 9×7 + 8 + 33 = **332bp**

Cutsite位置计算：
```python
target_start = 23 + 5  # 跳过Primer5和prefix
for i in range(10):
    t_start = target_start + i * (20 + 7)  # 每个target 20bp + PAM_Linker 7bp
    cutsite_start = t_start + 13           # 13bp conserved后
    cutsite_end = cutsite_start + 6        # 7bp cutsite
```

## 谱系示踪比对管线

完整处理流程：

```
  输入: (ref_seq, query_seq, cutsite列表)
         │
         ▼
  Step 1: 构建位置依赖gap惩罚数组
         │  build_gap_penalty_profile(ref_length, cutsites, ...)
         │  返回: (gap_open_profile, gap_extend_profile)
         │
         ▼
  Step 2: 位置感知DP比对
         │  affine_gap_alignment_position_aware(ref, query, profiles, ...)
         │  返回: (score, aligned_ref, aligned_query, stats)
         │
         ▼
  Step 3: 高密度mismatch→indel转换
         │  convert_dense_mismatch_to_indel(aligned_ref, aligned_query, ...)
         │  滑动窗口扫描，>34%密度+≥2错配→转换
         │
         ▼
  Step 4: 区域感知点突变过滤
         │  filter_point_mutations(aligned_ref, aligned_query, cutsites, ...)
         │  cutsite±3bp外矫正、gap邻域保护
         │
         ▼
  Step 5: 重新计分
         │  calculate_alignment_stats(final_ref, final_query)
         │  更新统计：矫正数、密集区转换标记
         │
         ▼
  输出: (score, final_ref, final_query, final_stats)
```

---

# 第三部分：跨靶点重复序列比对矫正

## 问题背景

CARLIN扩增子包含10个串联的靶点（Target），每个靶点由13bp保守区和7bp cutsite组成，靶点之间由7bp PAM-Linker（TGGAGTC）分隔。由于以下原因，扩增子内部存在大量重复序列：

1. **Target间共享序列**：多个Target包含完全相同的motif（如`ACAGTCG`出现在T1、T3、T9的cutsite区域）
2. **Linker重复**：9个PAM-Linker序列完全相同（`TGGAGTC`），造成整个扩增子的周期性结构
3. **短片段重复**：4-6bp的短序列（如`GACGA`、`ACTA`）在多个Target中重复出现

### DP算法的局限性

标准DP比对（含位置感知的谱系示踪算法）在处理重复序列时存在固有缺陷：

1. **随机匹配**：当参考序列中同一段序列出现在多个位置时，DP算法可能将query片段匹配到任意一个副本
2. **远端优先**：在仿射gap惩罚下，DP倾向于将短片段匹配到远端副本，因为远端匹配可避免开启新gap
3. **片段化**：本应连续的query序列被分割到两个远端副本，中间留下一段gap

典型问题模式：

```
参考序列: ...ACAGTCG... (T1正确位置) ...ACAGTCG... (T3远端副本)
query序列: ...ACAGTCG... (query实际对应T1) ...-------- (T3此处本应有连续序列)
                    ↓ DP算法错误匹配
query序列: ...-------- (T1正确位置空缺) ...ACAGTCG... (错误匹配到T3远端副本)
```

## 矫正算法

### 核心思想

利用生物学先验知识——**靶点区域的序列应该连续**，而非被分割到两个远端位置——在DP比对完成后，对比对结果进行后处理矫正。

### 算法流程

```
DP比对完成后的alignment
         │
         ▼
Step 1: 构建参考序列位置映射 (pos_map)
         │  将alignment列索引映射到参考序列位置
         │  pos_map[col] = ref_position; 插入列 = -1
         │
         ▼
Step 2: 扫描预定义的重复序列配置列表
         │  按序列长度从长到短处理（长序列优先）
         │
         ▼
Step 3: 对每个重复序列，检测错误匹配模式
         │  ① 远端副本有query碱基匹配
         │  ② 正确副本位置有gap空缺
         │  ③ 正确副本附近有query碱基（非真实删除）
         │
         ▼
Step 4: 搬迁碱基
         │  将query碱基从错误副本搬迁到正确副本
         │
         ▼
Step 5: 处理片段前导G碱基
         │  针对ACTGCACGACAGTCG的特殊后处理
         │
         ▼
输出矫正后的alignment
```

### 位置映射（pos_map）

pos_map数组将alignment的每一列映射到参考序列的对应位置：

```python
pos_map = []
cur = 0
for c in aligned_ref:
    pos_map.append(cur if c != '-' else -1)
    if c != '-':
        cur += 1
```

- 对于匹配/错配列（ar非`-`）：pos_map记录该列对应的参考序列位置
- 对于插入列（ar为`-`）：pos_map记录为`-1`

通过pos_map，可以精确找到重复序列各副本在alignment中的列索引，即使存在插入偏移。

### 重复序列配置

系统维护一个重复序列配置列表，每个条目包含重复序列本身和可选的特定(正确位置, 错误位置)对：

| 重复序列 | 长度 | 涉及靶点 | 特定位置对 |
|---------|------|---------|-----------|
| ACTGCACGACAGTCG | 15bp | T1↔T9 | 自动检测（首个副本为正确） |
| ACTCGCG | 7bp | T2↔T7 | 自动检测 |
| ACAGTCG | 7bp | T1↔T3↔T9 | (37,86), (37,253), (86,253) |
| GAGCGC | 6bp | T4↔T6 | 自动检测 |
| GCGACT | 6bp | T4↔T7 | 自动检测 |
| GATACG | 6bp | T5↔T10 | 自动检测 |
| ACGCAC | 6bp | T7↔T10 | 自动检测 |
| CGCGCA | 6bp | T2↔T5 | 自动检测 |
| CGACTA | 6bp | T4↔T9 | (109, 258) |
| GACGA | 5bp | T1↔T3 | 自动检测 |
| ACTA | 4bp | T9↔T3↔T4 | (260,83), (260,125), (125,83) |

**自动检测模式**（`specific_pairs=None`）：
- 使用`ref_seq.find(repeat)`查找所有副本
- 按参考序列上的出现顺序排列
- 以第一个副本为正确位置，其余为错误位置

**指定位置对模式**（`specific_pairs=[...]`）：
- 直接指定(正确位置, 错误位置)对
- 用于需要精确控制方向的情况（如ACTA从T9→T3而非T3→T9）

### 检测条件

每个位置对(correct_rp, wrong_rp)需要同时满足以下条件才执行搬迁：

1. **位置存在性**：correct_rp和wrong_rp处的rlen个碱基在alignment中都有对应的列（`pos_map`能定位到）

2. **参考序列匹配**：错误位置处的参考序列碱基与重复序列完全一致
   ```python
   all(ar[wrong_cols[k]] == repeat[k] for k in range(rlen))
   ```

3. **query匹配度**：错误位置处的query碱基与重复序列匹配度≥50%
   ```python
   match_cnt ≥ rlen * 0.5
   ```

4. **正确位置空缺**：正确位置处至少有50%的列为gap（未被其他搬迁填充）
   ```python
   gaps_at_correct ≥ rlen * 0.5
   ```

5. **邻域保护条件**：正确位置±3范围内有query碱基，避免搬入完全删除的区域
   ```python
   _adjacent_has_bases(correct_rp, rlen)
   ```

### 邻域保护函数

`_adjacent_has_bases(ref_pos, seg_len)`函数检查正确位置附近是否有query碱基，防止将碱基搬入一个完全被删除的靶点区域：

```python
def _adjacent_has_bases(ref_pos: int, seg_len: int) -> bool:
    for offset in range(ref_pos - 3, ref_pos + seg_len + 3):
        if 0 <= offset < len(ref_seq):
            col = find_column_by_position(pos_map, offset)
            if col is not None and aq[col] != '-':
                return True
    return False
```

检查范围是`[ref_pos-3, ref_pos+len+2]`（即正确位置前后各延伸3bp），确保目标区域不是完全删除。

### CGACTA的设计考量

`CGACTA`（6bp）在CARLIN参考序列中出现3次：
- ref[81:87]（T3起始附近，GACTACAGTCG中的前6bp）
- ref[109:115]（T4尾部）
- ref[258:264]（T9尾部，即ACTA前面2bp + ACTA 4bp）

使用自动检测模式会将ref[81]作为正确位置，导致本应在T9尾部的CGACTA（含ACTA）被错误地搬迁到T3起始处，破坏T9尾部的正确比对。

因此，CGACTA使用指定位置对`[(109, 258)]`，仅矫正T4↔T9之间的情况，避免影响T3起始区域。

### ACTA校正

`ACTA`（4bp）在CARLIN参考序列中出现3次：
- ref[83:87]（T3起始：GACTACAGTCG中的第2-5bp）
- ref[125:129]（T4尾部）
- ref[260:264]（T9末尾4bp）

DP算法可能将T9末尾的ACTA错误匹配到T3起始或T4尾部。指定位置对`[(260, 83), (260, 125), (125, 83)]`将从T3和T4向T9进行校正。

## 小片段跨Target矫正

### TAGTAT矫正

`TAGTAT`（6bp）是T8（ref[219:225]）的特有序列，但`GTCGAT`等相似序列可能出现在T9之后的Linker区域。DP算法可能将query中的TAGTAT匹配到T9附近区域。

矫正逻辑：
1. 检测T8的TAGTAT位置是否有≥4个gap（空缺）
2. 向后搜索50bp范围内是否有符合TAGTAT的query片段
3. 将query碱基搬迁回T8

### 单碱基A矫正

靶点T1倒数第4位（ref[44]）和T9倒数第4位（ref[260]）均为A碱基。DP算法可能将本应在T1的A匹配到T9附近。

矫正逻辑：
1. 检测T1的A位置是否为gap
2. 在T9的A位置±3范围内搜索孤立的A碱基
3. 将A搬迁回T1

## 孤立匹配清除

### 问题描述

当一个大片段的删除区域中间出现一个孤立的单碱基匹配时，该匹配会将该deletion片段分割成两个独立的缺失块。这在生物学上是不合理的——一个连续的删除事件不应被一个偶然匹配的碱基打断。

### 清除算法

遍历alignment，对每个位置i：

```python
# 检查条件
ar[i] != '-'            # 参考序列有碱基
aq[i] != '-'            # query有碱基
ar[i] == aq[i]          # 二者匹配（match）
ar[i-1] == '-' or aq[i-1] == '-'   # 左侧相邻列有gap
ar[i+1] == '-' or aq[i+1] == '-'   # 右侧相邻列有gap
# → 符合条件则将aq[i]转为'-'（即转为deletion）
```

该处理使两侧连续deletion合并，比对结果更简洁且更符合生物学实际。

## 完整矫正管线

谱系示踪比对的后处理管线依次执行：

```
DP比对原始输出
    │
    ▼
① correct_repetitive_misalignment  (跨靶点重复序列矫正)
    │  处理所有预定义的重复序列对
    │
    ▼
② correct_target_misalignments     (小片段跨靶点矫正)
    │  处理TAGTAT和单碱基A
    │
    ▼
③ remove_isolated_matches          (孤立匹配清除)
    │  清除打断连续deletion的孤立匹配
    │
    ▼
最终矫正结果
```

## 与标准算法的对比

| 特性 | 标准算法 | 谱系示踪算法 |
|------|---------|-------------|
| gap惩罚 | 全局固定 | 位置依赖 |
| 突变过滤 | 无 | cutsite区域感知 |
| indel检测 | 被动（DP决定） | 主动（密度检测+转换） |
| 假阳性控制 | 无 | 自动矫正保守区突变 |
| 结构先验 | 不使用 | 利用cutsite位置信息 |
| 适用场景 | 通用比对 | 谱系示踪/靶向基因编辑 |

## 性能分析

### 时间复杂度
- **标准算法**：$O(m \times n)$
- **谱系示踪算法**：$O(m \times n)$（增加的位置惩罚查询为$O(1)$，后处理为$O(L)$）

两者时间复杂度相同。

### 空间复杂度
- **标准算法**：$3 \times (m+1) \times (n+1)$ 个浮点数
- **谱系示踪算法**：额外增加 $2 \times m$ 个惩罚数组

### 执行时间估计
| 序列长度 | 标准算法 | 谱系示踪算法（含后处理） |
|----------|---------|----------------------|
| 100bp × 100bp | 5 ms | 6 ms |
| 332bp × 133bp | 45 ms | 50 ms |
| 500bp × 500bp | 125 ms | 140 ms |

谱系示踪算法增加约10-15%的运行时间，主要来自后处理步骤。

## 参考文献

1. Gotoh, O. (1982). An improved algorithm for matching biological sequences. *Journal of Molecular Biology*, 162(3), 705-708.
2. Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequences of two proteins. *Journal of Molecular Biology*, 48(3), 443-453.
3. Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.
