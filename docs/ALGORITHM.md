# 仿射gap惩罚序列比对算法原理

## 概述

本工具实现了基于Gotoh算法的仿射gap惩罚序列比对算法，是Needleman-Wunsch算法的仿射gap扩展。该算法专门为生物信息学序列比对设计，特别适用于CARLIN等基因编辑实验的序列分析。

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

## 默认参数设置

### CARLIN序列优化参数
| 参数 | 默认值 | 说明 | 设计意图 |
|------|--------|------|----------|
| `match_score` | 2.0 | 匹配得分 | 强化正确匹配 |
| `mismatch_penalty` | -3.0 | 错配惩罚 | 高于gap惩罚，减少错配 |
| `gap_open` | -2.0 | gap开启惩罚 | 相对较低，鼓励插入/缺失检测 |
| `gap_extend` | -0.1 | gap延伸惩罚 | 极低，鼓励长连续gap |

### 参数调整建议

#### 检测大片段缺失
```python
gap_open = -1.5      # 降低开启惩罚，鼓励检测缺失
gap_extend = -0.05   # 降低延伸惩罚，鼓励长连续gap
```

#### 严格匹配检测
```python
match_score = 3.0    # 提高匹配得分
mismatch_penalty = -5.0  # 增加错配惩罚
```

## 性能分析

### 时间复杂度
- **计算复杂度**：$O(m \times n)$
- **三个DP矩阵**：每个$(m+1) \times (n+1)$
- **回溯复杂度**：$O(m + n)$

### 空间复杂度
- **存储需求**：$3 \times (m+1) \times (n+1)$ 个浮点数
- **典型内存使用**：
  - 100bp × 100bp: ~0.8 MB
  - 500bp × 500bp: ~6 MB
  - 1000bp × 1000bp: ~24 MB
  - 332bp × 133bp (CARLIN): ~1.5 MB

### 执行时间估计
| 序列长度 | 近似时间 |
|----------|----------|
| 100bp × 100bp | 5 ms |
| 500bp × 500bp | 125 ms |
| 1000bp × 1000bp | 500 ms |
| 332bp × 133bp (CARLIN) | 45 ms |

## 算法验证

### 正确性验证
1. **边界条件测试**：完全匹配序列、完全不匹配序列
2. **gap检测测试**：不同长度和位置的插入/缺失
3. **参数敏感性测试**：验证参数调整对结果的影响

### 与MATLAB版本一致性
算法设计确保与MATLAB CARLIN版本的比对策略一致：
- 相同的gap惩罚模型
- 相同的回溯策略
- 相同的得分计算方式

## 实现细节

### 关键函数

#### `affine_gap_alignment()`
主比对函数，实现完整的DP算法和回溯。

#### `calculate_alignment_stats()`
计算比对统计信息：
- 匹配数、错配数
- gap数量和分布
- 相似度和一致性得分

#### `count_gap_blocks()`
统计连续gap区块，用于分析indel模式。

### 数值稳定性
- 使用`-np.inf`表示不可能状态
- 浮点数比较使用相对容差
- 避免数值下溢

## 参考文献

1. Gotoh, O. (1982). An improved algorithm for matching biological sequences. *Journal of Molecular Biology*, 162(3), 705-708.
2. Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443-453.
3. Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.

## 扩展建议

### 性能优化
1. **滚动数组**：将空间复杂度从$O(mn)$降低到$O(n)$
2. **并行计算**：利用多核CPU进行批量序列比对
3. **GPU加速**：使用CUDA实现大规模并行DP计算

### 功能扩展
1. **局部比对**：实现Smith-Waterman算法
2. **多序列比对**：扩展为渐进式多序列比对
3. **概率模型**：引入隐马尔可夫模型进行概率比对