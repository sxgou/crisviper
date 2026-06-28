# Algorithm / 算法原理

CrisViper 的核心是两个仿射 gap 比对模式（标准/谱系示踪），底层基于 Gotoh 算法，上层叠加 Numba JIT 编译加速的 DP 填充、NumPy 向量化预处理和后处理矫正管线。当 numba 不可用时自动降级为纯 NumPy/Python 实现。

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

**全局**（默认）：`M[0,0]=0`，其余首行首列 M 为 `-∞`，Ix 和 Iy 以 gap 惩罚初始化。回溯从右下角 `(m,n)` 走到左上角 `(0,0)`，剩余碱基填充 gap。两端 gap 按正常 gap 惩罚计分。不提供 --semi-global 选项。

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

## Backtrace / 回溯

DP 填充完成后，从右下角 `(m,n)` 回溯到 `(0,0)` 得到最优对齐路径。

### 终止状态选择

取三个状态矩阵在末位置的最大值作为终止得分：

```
max_score = max(M[m,n], Ix[m,n], Iy[m,n])
state = argmax([M[m,n], Ix[m,n], Iy[m,n]])
```

### 回溯规则

```
while i > 0 AND j > 0:
    if state == 'M':
        取 ref[i-1] 和 query[j-1] 对齐
        反查 M/Ix/Iy 中最大分数的前一状态（gap_exit_bonus 加入 Ix/Iy）
        i -= 1, j -= 1
    elif state == 'Ix':
        ref 侧填 gap '-'
        取 query[j-1] 对齐
        反查 M/Ix/Iy 中最大分数的前一状态（含 gap_open + gap_extend）
        j -= 1
    else:  # state == 'Iy'
        取 ref[i-1] 对齐
        query 侧填 gap '-'
        反查（含 gap_open + gap_extend）
        i -= 1

# 剩余列填充
while i > 0:  aligned_ref ← ref[i-1], aligned_query ← '-'
while j > 0:  aligned_ref ← '-', aligned_query ← query[j-1]
```

位置感知模式（`affine_gap_alignment_position_aware`）的回溯中，`gap_open_profile` 和 `gap_extend_profile` 根据 `(i,j)` 位置取 profile 值，使得 cutsite 区域的 gap 惩罚低于保守区，确保回溯路径与 DP 填充时的 position-aware 得分一致。

---

## Per-sequence pipeline / 单序列分析管线

`pipeline.py` 中 `align_single()` 按顺序执行以下步骤：

```
1. Full-length global alignment  (affine_gap_alignment / lineage_tracer_align)
2. Primer quality check          (_check_primer_quality)
3. Internal region extraction    (_extract_internal_region)
4. Background substitution corr. (correct_background_substitutions)
5. DEL→INS→DEL merge            (_merge_del_ins_del)
6. WT primer assembly            (_assemble_full_length)
7. Internal region re-scoring    (calculate_alignment_stats)
8. Mutation extraction           (extract_mutations)
9. Allele confidence filtering   (check_allele_confidence)
```

### 步骤说明

**1. 全局比对** — 标准模式调用 `affine_gap_alignment`，谱系模式调用 `lineage_tracer_align`（内部调用 `affine_gap_alignment_position_aware`）。详见前节 DP 递推。

**2. 引物质量检查** — 检查 5′ 和 3′ 引物区域的匹配碱基数是否达到阈值（默认 23 bp 中 ≥19，33 bp 中 ≥29）。引物锚定失败则标记为 anchor 失败。

**3. 内部区域提取** — 去除引物比对列，只保留靶标区域（internal region）用于下游突变分析。

**4. 背景点突变矫正** — 与旧版后处理不同，此矫正发生在 primer-trimmed 内部区域。详见下方 Background substitution correction 节。

**5. DEL→INS→DEL 合并** — 梯度 gap 惩罚在 cutsite 边界处可能将单个 INDEL 事件分裂为 DEL→INS→DEL 三段。此步骤检测这种碎片化模式并将碱基重新排列为连续 DEL + 连续 INS。

**6. WT 引物组装** — 用参考序列的野生型引物序列替换 query 的引物区域，避免引物区域的突变造成 allele 碎片化。

**7. 内部区域重新计分** — 在裁剪后的内部区域上重新计算比对统计量（match、mismatch、gap 数）。

**8. 突变提取** — 逐列扫描对齐序列，识别 substitution、deletion、insertion 和复合 INDEL 事件（`extract_mutations`）。相邻的插入+删除被贪婪合并为 INDEL 复合事件。

**9. Allele 置信度过滤** — 纯点突变 allele 需 read 数 ≥ `--min-reads-sub`（默认 5），含 indel 的 allele 需 ≥ `--min-reads-indel`（默认 0）。野生型直接通过。

### Background substitution correction / 背景点突变矫正

高通量测序中的 PCR 错误和测序错误会在非编辑位点引入假阳性点突变。背景矫正算法自动识别并过滤这些假阳性：

1. 对每条序列，遍历其比对列中的点突变
2. 排除 cutsite 窗口（`--sub-window`，默认 ±3 bp）内的突变
3. 排除紧邻 indel 的突变（`--keep-sub-indel-window`，默认 ±3 bp）——这些突变可能参与了修复过程
4. 剩余点突变视为测序/PCR 噪声，矫正回 reference 碱基

可通过 `--correct-bg-sub` 关闭此功能。

### 上游 DEPRECATED 算法说明（MATLAB 遗留）

以下后处理步骤存在于 MATLAB 原始实现（`reports/generate_text_output.m`）但已被 Python 实现替换为 DP 原生特征或 pipeline 内步骤，仅供参考：

- **跨靶点重复序列矫正**: 旧版 MATLAB 的后处理步骤，处理 CARLIN 重复 motif 导致的错配。Python 版通过 DP 原生特征（homology_penalty + isolated_base_penalty）在比对阶段解决。
- **孤立匹配清除**: 旧版后处理。Python 版通过 `isolated_base_penalty`（DP 原生特征）在比对阶段抑制。
- **密集错配 → indel 转换**: 旧版后处理（滑动窗口替换）。Python 版通过 `dense_mismatch_penalty` 在 DP 填充阶段（`_build_score_matrix` + `affine_gap_alignment_position_aware`）直接修改替换得分矩阵，**不是**独立后处理步骤。
- **点突变过滤（cutsite 窗口外）**: 此功能已整合进 Background substitution correction。

### Event-level aggregation / 事件级聚合

后处理完成后，结果聚合成事件级统计表（`event_level_details.tsv`）：每个突变事件由 `(type, start_pos, length)` 唯一标识，统计序列数、reads 数、覆盖的 target（27bp 窗口：保守区 13bp + cutsite + linker 7bp）、跨越 target 数量和起止 target 名称。

---

## Numba JIT 编译加速

核心 DP 填充函数通过 Numba `@jit(nopython=True)` 编译为机器码，提供约 **28×** 加速（v1.1.0）。

### JIT 函数列表

| 函数 | 位置 | 用途 |
|------|------|------|
| `_dp_fill_numba` | alignment.py:197-239 | 位置感知模式 DP 填充（谱系模式） |
| `_dp_fill_numba_standard` | alignment.py:242-291 | 标准模式 DP 填充 |
| `_compute_run_len_numba` | alignment.py:136-142 | match 连续长度反向扫描 |
| `_dense_mm_density_numba` | alignment.py:151-187 | 密集错配密度对角线滑动窗口 |

### 降级机制

当 numba 不可用时（`ImportError`），自动降级为纯 Python 实现：

```python
try:
    from numba import jit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    def jit(nopython=True, cache=True):
        def wrapper(func):
            return func
        return wrapper
```

每个 JIT 函数通过一个 Python wrapper 间接调用（如 `_compute_run_len` → `_compute_run_len_numba`），便于测试和 mock。

### 编译缓存

Numba 将编译后的函数缓存为 `.nbc`/`.nbi` 文件，后续运行跳过重编译。

## Vectorization / 向量化与预处理

### NumPy 向量化预处理

DP 填充之前，`_build_score_matrix` 用 NumPy 向量化操作完成以下预处理：

1. **替换得分矩阵**: 将 ref/query 编码为 int8 数组，通过广播比较 `ref_codes[:, None] == qry_codes[None, :]` 一次生成 m×n 的 match/mismatch 矩阵
2. **同源惩罚**: 对检测到的同源区域，通过 `sub_score + hp[:, None] * is_match` 向量化减去惩罚
3. **短匹配折扣**: `_compute_run_len_numba`（Numba JIT）反向扫描计算 match 连续长度，短于阈值的区域 `sub_score[mask] = match_score * short_match_discount`
4. **孤立碱基惩罚**: `iso_mask = run_len[:m, :n] == 1` 定位孤立 match，`sub_score[iso_mask] += isolated_base_penalty`
5. **密集错配惩罚**: `_dense_mm_density_numba`（Numba JIT）沿对角线滑动窗口计算密度，超出阈值时 `sub_score[dense_mask] += dense_mismatch_penalty`

### Numba JIT DP 填充（主路径）

DP 填充矩阵 `M`、`Ix`、`Iy` 的逐行计算由 Numba JIT 编译执行（纯 Python 嵌套循环，Numba 编译为 LLVM IR → 机器码）：

```python
# _dp_fill_numba 核心循环（简写）
for i in range(1, m + 1):
    for j in range(1, n + 1):
        Iy[i,j] = max(M[i-1,j]+go+ge, Ix[i-1,j]+go+ge, Iy[i-1,j]+ge)
        best = max(M[i-1,j-1], Ix[i-1,j-1]+geb, Iy[i-1,j-1]+geb)
        M[i,j] = sub_score[i-1,j-1] + best
        Ix[i,j] = max(M[i,j-1]+go+ge, Ix[i,j-1]+ge, Iy[i,j-1]+go+ge)
```

### NumPy 向量化降级路径（Numba 不可用时）

当 numba 不可用时，使用纯 NumPy 向量化的降级路径。Ix 在同列中存在数据依赖，无法向量化；Iy 和 M 只依赖前一行的值，整行可用一次 NumPy 操作完成：

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

### Performance / 性能

**Numba JIT 加速效果（v1.1.0）：** 单序列 4.8 ms（Numba JIT）vs 136 ms（纯 NumPy 降级），加速比约 **28×**。

**历史性能对比：**

| 版本 | 单序列 | 500条×12线程 | 23430条×12线程 |
|------|--------|-------------|----------------|
| v1.0.0（纯 Python） | 1.60 s | 34 s | ~38 min |
| v1.0.x（NumPy 向量化） | 0.37 s | 19 s | ~15 min |
| v1.1.0+（Numba JIT） | **4.8 ms** | — | — |

批处理加速比低于单序列，原因是进程间通信开销和负载不均。以上性能数据在 24 核测试机器上测得。线程数建议不超过 CPU 物理核心数（推荐 max 12，超出后因进程间通信和内存带宽限制收益递减）。实际加速效果因 CPU、序列长度和数量而异。

---

## Downstream analysis modules / 下游分析模块

比对管线完成后，以下模块可进一步分析结果：

### Denoiser（crisviper/denoiser.py）

UMI/CB（Unique Molecular Identifier / Cell Barcode）去噪。使用 directional adjacency clustering 对相近序列进行去重，降低 PCR 重复和测序噪声。

### Caller（crisviper/caller.py）

Allele calling（等位基因判定）。提供两个策略：
- **coarse-grain**: 按突变指纹聚类，合并相同编辑模式的 allele
- **exact**: 精确匹配，仅合并序列完全一致的 allele

### Metrics（crisviper/metrics.py）

多样性/异质性度量。计算 Shannon entropy、effective alleles 等群体遗传学指标，评估编辑结果的克隆复杂度。

### Threshold（crisviper/threshold.py）

UMI/CB 过滤阈值计算。基于 read count 分布自动计算合适的过滤阈值。

这些模块可通过 YAML 配置文件（`--config`）中的 `pipeline` 段启用和配置，也可作为编程接口直接调用。

---

## References / 参考文献

1. Gotoh, O. (1982). An improved algorithm for matching biological sequences. *J. Mol. Biol.*, 162(3), 705–708.
2. Needleman, S. B. & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequences of two proteins. *J. Mol. Biol.*, 48(3), 443–453.
3. Smith, T. F. & Waterman, M. S. (1981). Identification of common molecular subsequences. *J. Mol. Biol.*, 147(1), 195–197.
