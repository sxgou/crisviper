# Algorithm / 算法原理

CrisViper 的核心是两个仿射 gap 比对模式（标准/谱系示踪），底层基于 Gotoh 算法，上层叠加 NumPy 向量化 DP 和后处理矫正管线。

---

## Standard alignment / 标准比对

### 仿射 gap 模型

Gap 成本不是线性的，而是由开启（open）和延伸（extend）两部分组成：总成本 = $g_o + k \cdot g_e$。这比线性 gap 更符合生物学实际——一个连续的长 indel 比多个分散的短 indel 更常见。

### DP 递推

三个状态矩阵，序列 $X$（长度 $m$）和 $Y$（长度 $n$）：

**Match/mismatch** $M(i,j)$:

$$
M(i,j) = s(x_i, y_j) + \max\begin{cases}
M(i-1, j-1) \\
I_x(i-1, j-1) \\
I_y(i-1, j-1)
\end{cases}
$$

**Gap in X (insertion in Y)** $I_x(i,j)$:

$$
I_x(i,j) = \max\begin{cases}
M(i, j-1) + g_o + g_e \\
I_x(i, j-1) + g_e \\
I_y(i, j-1) + g_o + g_e
\end{cases}
$$

**Gap in Y (deletion in X)** $I_y(i,j)$:

$$
I_y(i,j) = \max\begin{cases}
M(i-1, j) + g_o + g_e \\
I_x(i-1, j) + g_o + g_e \\
I_y(i-1, j) + g_e
\end{cases}
$$

$s(x_i, y_j)$ 是替换得分矩阵，match_score 或 mismatch_penalty。

**半全局**（默认）：首行首列初始化为 0，不惩罚序列两端的 gap。仅支持半全局模式，不提供 --global 选项。

---

## Lineage mode / 谱系示踪模式

### 设计动机

多靶点谱系示踪扩增子（如 CARLIN，332 bp，10 个串联 target）有明确的生物学先验：indel 几乎只发生在 cutsite，不会发生在保守区。标准 Gotoh 算法用固定的全局 gap 惩罚，无法利用这个先验。谱系模式把它嵌入 DP 本身。

### Gradient-based gap penalties / 梯度 gap 惩罚

gap_open 和 gap_extend 乘以位置因子 $s(p)$ 来调整：

$$
g_o(p) = g_o^{\text{base}} \times s(p), \quad g_e(p) = g_e^{\text{base}} \times s(p)
$$

与传统的离散区域划分不同，谱系模式使用基于 smoothstep 的**连续梯度**。每个 cutsite 的位置因子由三部分叠加决定：

1. **Cutsite 中心**：惩罚最低（`--min-scale`，默认 1.0），鼓励在切割位点开 gap
2. **边界衰减**：以平滑步进函数 $\operatorname{smoothstep}(t)=3t^2-2t^3$ 在 $[0,1]$ 区间内从最小惩罚过渡到最大惩罚，渐变半径由 `--gradient-radius` 控制（默认自动计算）
3. **Cutsite 边界**：叠加 `--cutsite-edge-scale`（默认 2.0）倍率的平滑峰，抑制边界附近的不精确 indel

保守区惩罚平滑过渡到 `--max-scale`（默认 6.0），避免离散区域边界处的突变效应。

| 位置 | $s(p)$ | 效果 |
|------|--------|------|
| Cutsite 中心 | 1.0 | 鼓励开 gap |
| Cutsite 边界 (±σ) | ~2.0 | 中等抑制 |
| 保守区（远离 cutsite） | 6.0 | 强烈抑制 |

### Cut site detection / 切割位点检测

标准 CARLIN 结构（332 bp，10 target）会自动推断：每个 target = 13 bp 保守区 + 7 bp cutsite，之间由 7 bp PAM-Linker 分隔。

非标准结构通过 `--cutsites` 传入 JSON：

```json
[
  {"name": "Target1", "start": 41, "end": 47}
]
```

坐标 0-based、inclusive。注：也接受 `{"cutsites": [...]}` 包装格式。

### DP-native features / DP 原生特征

这些特征直接修改 DP 递推式，让比对器在比对过程中做出更好的决策，而不是依赖后处理来修正。

**Gap exit strength.** 从 gap 状态切回 match 时额外施加惩罚。效果是抑制碎片化的短 match 区域，将小的 match 片段"推入"相邻 gap，合并为连续 indel。与 `isolated_base_penalty` 协同：前者惩罚 gap→M 的过渡，后者惩罚过渡后只有 1bp 匹配的场景。数学上，在 M 递推的 max 中，来自 Ix/Iy 的分支加上 `gap_exit_strength`（负值）：

$$
M(i,j) = s(x_i, y_j) + \max\begin{cases}
M(i-1, j-1) \\
I_x(i-1, j-1) + \text{gap\_exit\_bonus} \\
I_y(i-1, j-1) + \text{gap\_exit\_bonus}
\end{cases}
$$

**Short match discount.** 预处理阶段用反向扫描为每个单元格计算 match 连续长度。短于 `--short-match-window` 的区域，match_score 被 `--short-match-discount` 系数打折。这防止了重复序列区域中短片段偶然匹配锚定比对结果。

**Dense mismatch penalty.** 沿 DP 矩阵对角线用 cumsum 计算滑动窗口内的 mismatch 密度（见向量化章节）。密度超出阈值时，对替换得分施加额外惩罚，使 DP 倾向开 gap 而不是累积错配。

**Homology penalty.** 扫描参考序列的局部自相似性，在同源区域降低 match_score。抑制 DP 将 reads 比对到重复元件的错误副本。

**Isolated base penalty.** 单碱基匹配两侧都是 gap，会把一个连续 deletion 分成两个不连续的片段——生物学上不太可能。孤立碱基惩罚降低这类 match 的得分，让 DP 倾向选择更长的连续 deletion。

---

## Post-alignment correction / 后处理矫正

按顺序执行：

```
DP output
  → 跨靶点重复序列矫正 (11 组, 4–15 bp)
  → 小片段跨靶点矫正 (TAGTAT, 单碱基 A)
  → 孤立匹配清除
  → 密集错配 → indel 转换 (后处理备选)
  → 点突变过滤 (cutsite 窗口外)
  → 背景点突变矫正 (indel 邻近窗口保留)
  → 事件级聚合 (event_level_details.tsv)
```

### Repetitive misalignment correction

CARLIN 的 10 个串联 target 中包含大量重复 motif（如 `ACAGTCG` 同时出现在 T1、T3、T9；9 个完全相同的 `TGGAGTC` PAM-Linker）。DP 可能把 query 片段匹配到错误副本，产生分裂的 alignment。矫正算法检测这种模式——错误副本处有 query 碱基匹配、正确副本处为 gap——并将碱基搬回原位。

处理 11 组重复序列，按长度从长到短：`ACTGCACGACAGTCG` (15 bp)、`ACTCGCG`/`ACAGTCG` (7 bp)、6 个 6 bp 基序、`GACGA` (5 bp)、`ACTA` (4 bp)。

每个位置对需同时满足：错误位置处 reference 和 query 与重复序列匹配度 ≥ 50%、正确位置处 ≥ 50% 的列为 gap、正确位置邻域有 query 碱基（防止搬入完全删除的区域）。

### Isolated match removal / 孤立匹配清除

被 gap 夹在中间的单碱基 match 把一个 deletion 切成两块。算法将其转为 gap，合并 deletion。

### Dense mismatch → indel conversion

滑动窗口（默认 6 bp）内 mismatch 密度 > 34% 且至少 2 个错配时，窗口内的 reference 碱基替换为 gap（将对应区域视为 query 的 insertion）。

### Point mutation filtering / 点突变过滤

cutsite ±`sub_window` bp 范围外的点突变被矫正回 reference 碱基。紧邻 gap 的突变不受此规则影响——这些突变可能参与了修复过程，保留它们避免假阴性。

### Background substitution correction / 背景点突变矫正

高通量测序中的 PCR 错误和测序错误会在非编辑位点引入假阳性点突变。背景矫正算法自动识别并过滤这些假阳性：

1. 对每条序列，遍历其所有点突变
2. 排除 cutsite 窗口（`--sub-window`，默认 ±3 bp）内的突变
3. 排除紧邻 indel 的突变（`--keep-sub-indel-window`，默认 ±3 bp）——这些突变可能参与了修复过程
4. 剩余点突变视为测序/PCR 噪声，矫正回 reference 碱基

可通过 `--correct-bg-sub` 关闭此功能。

### Event-level aggregation / 事件级聚合

后处理完成后，结果聚合成事件级统计表（`event_level_details.tsv`）：每个突变事件由 `(type, start_pos, length)` 唯一标识，统计序列数、reads 数、覆盖的 target（20bp 窗口：保守区 13bp + cutsite 7bp）、跨越 target 数量和起止 target 名称。

---

## Vectorization / 向量化

Ix 在同列中存在数据依赖（依赖于 Ix 前一列的值），无法向量化。但 Iy 和 M 只依赖前一行的值，整行可以用一次 NumPy 操作完成：

```python
prev_best = np.maximum(M[i-1], np.maximum(Ix[i-1], Iy[i-1]))
M[i, 1:] = score_row + prev_best[:n-1]

ge_arr = gap_extend_profile[i]
go_ge_arr = gap_open_profile[i] + ge_arr
Iy[i, 1:] = np.maximum(np.maximum(M[i-1, 1:] + go_ge_arr,
                                   Ix[i-1, 1:] + go_ge_arr),
                        Iy[i-1, 1:] + ge_arr)
```

由于 $n \leq 332$（参考序列长度），Ix 的 Python 顺序循环开销可以忽略。

**Density via cumsum.** 密集错配密度计算从 $O(m \times n \times w)$ 的三重循环优化为纯 NumPy 的 $O(m \times n)$。方法：构建对角线索引矩阵 → 沿对角线排序 → cumsum → 滑动窗口差值。

### Performance / 性能

| | 优化前 | 优化后 | 加速比 |
|--|--------|-------|--------|
| 单序列（全部特性）| 1.60 s | 0.37 s | 4.3× |
| 500 条（12 线程）| 34 s | 19 s | 1.8× |
| 23430 条（12 线程）| ~38 min | ~15 min | 2.5× |

批处理加速比低于单序列，原因是进程间通信开销和负载不均。以上性能数据在 24 核测试机器上使用 `--threads 12` 测得，软件本身没有线程数上限。实际加速效果因 CPU、序列长度和数量而异。

---

## References / 参考文献

1. Gotoh, O. (1982). An improved algorithm for matching biological sequences. *J. Mol. Biol.*, 162(3), 705–708.
2. Needleman, S. B. & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequences of two proteins. *J. Mol. Biol.*, 48(3), 443–453.
3. Smith, T. F. & Waterman, M. S. (1981). Identification of common molecular subsequences. *J. Mol. Biol.*, 147(1), 195–197.
