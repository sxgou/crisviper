#!/usr/bin/env python3
"""
仿射gap惩罚的精确序列比对算法
基于Gotoh算法（Needleman-Wunsch的仿射gap扩展）

算法特点：
1. 最小化错配（错配惩罚高于gap惩罚）
2. 减少短gap（gap extension惩罚极低，鼓励连续长gap）
3. 鼓励插入而非错配（gap opening惩罚低于错配惩罚）
4. 支持半全局比对（两端gap不惩罚）

作者：Hermes Agent
日期：2026-04-23
"""

import numpy as np
from typing import Tuple, Dict, List, Optional
from dataclasses import dataclass, field


def affine_gap_alignment(ref_seq: str, query_seq: str,
                         match_score: float = 2.0,
                         mismatch_penalty: float = -3.0,
                         gap_open: float = -2.0,
                         gap_extend: float = -0.1,
                         semi_global: bool = True,
                         fit_mode: bool = False) -> Tuple[float, str, str, Dict]:
    """
    仿射gap惩罚的精确序列比对算法

    参数：
        ref_seq: reference序列（字符串）
        query_seq: query序列（字符串）
        match_score: 匹配得分（正数）
        mismatch_penalty: 错配惩罚（负数）
        gap_open: gap开启惩罚（负数）
        gap_extend: gap延伸惩罚（负数，绝对值小于gap_open）
        semi_global: 是否使用半全局比对（两端gap不惩罚）
        fit_mode: 短序列fit模式（query必须完整使用，ref两端自由）
                  优先级高于semi_global

    返回：
        score: 比对得分
        aligned_ref: 比对后的reference序列（含gap）
        aligned_query: 比对后的query序列（含gap）
        stats: 统计信息字典
    """
    m, n = len(ref_seq), len(query_seq)

    # 三个DP矩阵
    M = np.zeros((m + 1, n + 1), dtype=float)  # 匹配/错配状态
    Ix = np.zeros((m + 1, n + 1), dtype=float)  # gap在ref中（query插入）
    Iy = np.zeros((m + 1, n + 1), dtype=float)  # gap在query中（ref缺失）

    # 初始化
    if fit_mode:
        # fit模式：query必须完整使用，ref两端自由
        # ref可以在任意位置开始和结束对应query，但query不能头尾缺失
        M[0, 0] = 0.0
        Ix[0, 0] = -np.inf
        Iy[0, 0] = -np.inf

        for i in range(1, m + 1):
            M[i, 0] = 0.0           # ref可以自由开始
            Ix[i, 0] = -np.inf
            Iy[i, 0] = gap_open + gap_extend * (i - 1)  # ref缺失（有代价）

        for j in range(1, n + 1):
            M[0, j] = -np.inf       # query不能在ref之前开始
            Ix[0, j] = -np.inf
            Iy[0, j] = -np.inf
    elif semi_global:
        # 半全局：第一行和第一列初始化为0（两端gap不惩罚）
        M[0, 0] = 0.0
        Ix[0, 0] = -np.inf
        Iy[0, 0] = -np.inf
        
        for i in range(1, m + 1):
            M[i, 0] = 0.0
            Ix[i, 0] = -np.inf
            Iy[i, 0] = 0.0  # 允许ref缺失而不惩罚
            
        for j in range(1, n + 1):
            M[0, j] = 0.0
            Ix[0, j] = 0.0  # 允许query插入而不惩罚
            Iy[0, j] = -np.inf
    else:
        # 全局比对
        M[0, 0] = 0.0
        Ix[0, 0] = gap_open
        Iy[0, 0] = gap_open
        
        for i in range(1, m + 1):
            M[i, 0] = -np.inf
            Ix[i, 0] = -np.inf
            Iy[i, 0] = gap_open + gap_extend * (i - 1)
            
        for j in range(1, n + 1):
            M[0, j] = -np.inf
            Ix[0, j] = gap_open + gap_extend * (j - 1)
            Iy[0, j] = -np.inf
    
    # 填充DP矩阵（使用行预取和Python内置max加速）
    go_ge = gap_open + gap_extend
    for i in range(1, m + 1):
        Mi_1 = M[i - 1]; Mi = M[i]
        Ixi_1 = Ix[i - 1]; Ixi = Ix[i]
        Iyi_1 = Iy[i - 1]; Iyi = Iy[i]
        r_char = ref_seq[i - 1]
        for j in range(1, n + 1):
            s = match_score if r_char == query_seq[j - 1] else mismatch_penalty
            # M状态
            a, b, c = Mi_1[j - 1], Ixi_1[j - 1], Iyi_1[j - 1]
            Mi[j] = s + (a if a >= b and a >= c else (b if b >= c else c))
            # Ix状态
            a, b, c = Mi[j - 1] + go_ge, Ixi[j - 1] + gap_extend, Iyi[j - 1] + go_ge
            Ixi[j] = a if a >= b and a >= c else (b if b >= c else c)
            # Iy状态
            a, b, c = Mi_1[j] + go_ge, Ixi_1[j] + go_ge, Iyi_1[j] + gap_extend
            Iyi[j] = a if a >= b and a >= c else (b if b >= c else c)
    
    # 确定最终得分和回溯起点
    if fit_mode:
        # fit模式：从最后一列(j=n)中取最大值（query必须完全使用）
        max_score = -np.inf
        max_i, max_j, max_state = 0, n, 'M'
        for i in range(m + 1):
            for state, val in [('M', M[i, n]), ('Ix', Ix[i, n]), ('Iy', Iy[i, n])]:
                if val > max_score:
                    max_score = val
                    max_i, max_j, max_state = i, n, state
    elif semi_global:
        # 半全局：取最后一行和最后一列的最大值（两端gap不惩罚）
        max_score = -np.inf
        max_i, max_j, max_state = 0, 0, 'M'

        for i in range(m + 1):
            for state, val in [('M', M[i, n]), ('Ix', Ix[i, n]), ('Iy', Iy[i, n])]:
                if val > max_score:
                    max_score = val
                    max_i, max_j, max_state = i, n, state

        for j in range(n + 1):
            for state, val in [('M', M[m, j]), ('Ix', Ix[m, j]), ('Iy', Iy[m, j])]:
                if val > max_score:
                    max_score = val
                    max_i, max_j, max_state = m, j, state
    else:
        # 全局比对：右下角
        max_score = max(M[m, n], Ix[m, n], Iy[m, n])
        max_i, max_j = m, n
        max_state = 'M' if M[m, n] >= max(Ix[m, n], Iy[m, n]) else \
                   ('Ix' if Ix[m, n] >= Iy[m, n] else 'Iy')

    # 回溯
    i, j = max_i, max_j
    state = max_state
    aligned_ref, aligned_query = [], []

    while (i > 0 or j > 0) and (fit_mode or semi_global or (i > 0 and j > 0)):
        if state == 'M':
            if i > 0 and j > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append(query_seq[j - 1])
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
                j -= 1
            elif i > 0:
                # ref延伸超出query → 作为deletion
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                i -= 1
                state = 'Iy'
            else:
                break
        elif state == 'Ix':
            aligned_ref.append('-')
            aligned_query.append(query_seq[j - 1])
            # 前一个状态
            scores = [
                M[i, j - 1] + gap_open + gap_extend,
                Ix[i, j - 1] + gap_extend,
                Iy[i, j - 1] + gap_open + gap_extend
            ]
            prev = np.argmax(scores)
            state = ['M', 'Ix', 'Iy'][prev]
            j -= 1
        else:  # state == 'Iy'
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append('-')
            # 前一个状态
            scores = [
                M[i - 1, j] + gap_open + gap_extend,
                Ix[i - 1, j] + gap_open + gap_extend,
                Iy[i - 1, j] + gap_extend
            ]
            prev = np.argmax(scores)
            state = ['M', 'Ix', 'Iy'][prev]
            i -= 1
    
    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))
    
    # 统计信息
    stats = calculate_alignment_stats(aligned_ref, aligned_query)
    stats['score'] = max_score
    
    return max_score, aligned_ref, aligned_query, stats


def calculate_alignment_stats(aligned_ref: str, aligned_query: str) -> Dict:
    """
    计算比对统计信息
    
    参数：
        aligned_ref: 比对后的reference序列
        aligned_query: 比对后的query序列
    
    返回：
        统计信息字典
    """
    matches = sum(1 for a, b in zip(aligned_ref, aligned_query) if a == b and a != '-')
    mismatches = sum(1 for a, b in zip(aligned_ref, aligned_query) 
                     if a != b and a != '-' and b != '-')
    gaps_in_ref = aligned_ref.count('-')
    gaps_in_query = aligned_query.count('-')
    
    # 计算连续gap统计
    gap_blocks_ref = count_gap_blocks(aligned_ref)
    gap_blocks_query = count_gap_blocks(aligned_query)
    
    # 计算平均gap长度
    avg_gap_len_ref = np.mean(gap_blocks_ref) if gap_blocks_ref else 0
    avg_gap_len_query = np.mean(gap_blocks_query) if gap_blocks_query else 0
    
    # 计算相似度
    alignment_length = len(aligned_ref)
    similarity = matches / alignment_length if alignment_length > 0 else 0
    
    stats = {
        'matches': matches,
        'mismatches': mismatches,
        'gaps_in_ref': gaps_in_ref,
        'gaps_in_query': gaps_in_query,
        'gap_blocks_ref': gap_blocks_ref,
        'gap_blocks_query': gap_blocks_query,
        'avg_gap_len_ref': avg_gap_len_ref,
        'avg_gap_len_query': avg_gap_len_query,
        'alignment_length': alignment_length,
        'similarity': similarity,
        'identity': matches / (matches + mismatches) if (matches + mismatches) > 0 else 0
    }
    
    return stats


def count_gap_blocks(seq: str) -> List[int]:
    """
    统计序列中的连续gap区块长度
    
    参数：
        seq: 比对序列（可能包含gap字符'-'）
    
    返回：
        连续gap区块长度列表
    """
    blocks = []
    in_gap = False
    current_length = 0
    
    for c in seq:
        if c == '-':
            if not in_gap:
                in_gap = True
                current_length = 1
            else:
                current_length += 1
        else:
            if in_gap:
                blocks.append(current_length)
                in_gap = False
    
    if in_gap:
        blocks.append(current_length)
    
    return blocks


def print_alignment(aligned_ref: str, aligned_query: str, width: int = 100):
    """
    打印比对结果
    
    参数：
        aligned_ref: 比对后的reference序列
        aligned_query: 比对后的query序列
        width: 每行显示的字符数
    """
    print("\n比对结果:")
    for i in range(0, len(aligned_ref), width):
        ref_slice = aligned_ref[i:i + width]
        query_slice = aligned_query[i:i + width]
        
        # 添加行号
        print(f"Ref  {i:4d}: {ref_slice}")
        print(f"Query{i:4d}: {query_slice}")
        
        # 添加匹配指示行
        match_line = " " * 9
        for r, q in zip(ref_slice, query_slice):
            if r == q and r != '-':
                match_line += '|'
            elif r != '-' and q != '-':
                match_line += '.'
            else:
                match_line += ' '
        print(match_line)
        print()


def print_stats(stats: Dict):
    """
    打印统计信息
    
    参数：
        stats: 统计信息字典
    """
    print("\n比对统计:")
    print(f"  比对得分: {stats['score']:.2f}")
    print(f"  比对长度: {stats['alignment_length']}")
    print(f"  匹配数: {stats['matches']}")
    print(f"  错配数: {stats['mismatches']}")
    print(f"  相似度: {stats['similarity']:.2%}")
    print(f"  一致性: {stats['identity']:.2%}")
    print(f"  Reference中gap数: {stats['gaps_in_ref']}")
    print(f"  Query中gap数: {stats['gaps_in_query']}")
    
    if stats['gap_blocks_ref']:
        print(f"  Reference中gap区块: {stats['gap_blocks_ref']}")
        print(f"  Reference平均gap长度: {stats['avg_gap_len_ref']:.2f}")
    
    if stats['gap_blocks_query']:
        print(f"  Query中gap区块: {stats['gap_blocks_query']}")
        print(f"  Query平均gap长度: {stats['avg_gap_len_query']:.2f}")


# ═══════════════════════════════════════════════════════════════
# 结构感知比对 — 用于谱系示踪实验的靶向分析
# ═══════════════════════════════════════════════════════════════

@dataclass
class CutsiteRegion:
    """cutsite区域定义 (0-indexed, inclusive)"""
    name: str               # 如 "Target1"
    start: int              # cutsite起始位置 (参考序列上的坐标)
    end: int                # cutsite结束位置 (参考序列上的坐标)

@dataclass
class StructureConfig:
    """扩增子结构配置"""
    cutsites: List[CutsiteRegion] = field(default_factory=list)
    primers: Dict[str, Tuple[int, int]] = field(default_factory=dict)


def build_gap_penalty_profile(
    ref_length: int,
    cutsites: List[CutsiteRegion],
    base_gap_open: float = -2.0,
    base_gap_extend: float = -0.1,
    cutsite_scale: float = 1.0,
    flank_scale: float = 2.0,
    far_scale: float = 6.0,
    flank_width: int = 3
) -> Tuple[np.ndarray, np.ndarray]:
    """
    构建位置依赖的gap惩罚数组。

    设计原则:
      cutsite区域: 基准惩罚 (1x) — 鼓励gap开启
      cutsite周边 (flanking): 2x惩罚 — 中等抑制
      保守/骨架区域: far_scale 惩罚 — 强烈抑制gap

    参数:
        ref_length: 参考序列长度
        cutsites: cutsite区域列表
        base_gap_open: 基准gap开启惩罚 (负值)
        base_gap_extend: 基准gap延伸惩罚 (负值)
        cutsite_scale: cutsite内惩罚倍率 (默认1.0)
        flank_scale: 侧翼区惩罚倍率 (默认2.0)
        far_scale: 远离cutsite区域惩罚倍率 (默认6.0)
        flank_width: 侧翼区宽度(bp) (默认3)

    返回:
        (gap_open_profile, gap_extend_profile): 长度为ref_length的惩罚数组
    """
    # 初始化: 所有位置使用 far_scale
    gap_open_profile = np.full(ref_length, base_gap_open * far_scale)
    gap_extend_profile = np.full(ref_length, base_gap_extend * far_scale)

    # 标记每个位置的最低有效倍率
    effective_scale = np.full(ref_length, far_scale)

    for cs in cutsites:
        cs_start = max(0, cs.start)
        cs_end = min(ref_length - 1, cs.end)

        # cutsite内部: 基准倍率
        for pos in range(cs_start, cs_end + 1):
            effective_scale[pos] = min(effective_scale[pos], cutsite_scale)

        # 侧翼区域 (cutsite两侧 flank_width bp): flank_scale
        for offset in range(1, flank_width + 1):
            left_pos = cs_start - offset
            if left_pos >= 0:
                effective_scale[left_pos] = min(effective_scale[left_pos], flank_scale)
            right_pos = cs_end + offset
            if right_pos < ref_length:
                effective_scale[right_pos] = min(effective_scale[right_pos], flank_scale)

    # 根据有效倍率计算惩罚值
    for pos in range(ref_length):
        gap_open_profile[pos] = base_gap_open * effective_scale[pos]
        gap_extend_profile[pos] = base_gap_extend * effective_scale[pos]

    return gap_open_profile, gap_extend_profile


def affine_gap_alignment_position_aware(
    ref_seq: str,
    query_seq: str,
    gap_open_profile: np.ndarray,
    gap_extend_profile: np.ndarray,
    match_score: float = 2.0,
    mismatch_penalty: float = -3.0,
    semi_global: bool = True
) -> Tuple[float, str, str, Dict]:
    """
    位置感知的仿射gap序列比对。

    使用位置依赖的gap惩罚数组，允许在cutsite区域优先开启gap。
    gap_open_profile[i] 和 gap_extend_profile[i] 对应参考序列第i位的惩罚。

    semi_global模式:
      True  (默认) — 半全局，首尾ref gap免罚，query可在任意位置起止
      False — 全局，双端严格，seq必须耗尽
      'anchor5' — 锚定5'端，query必须从ref[0]开始

    返回: (score, aligned_ref, aligned_query, stats)
    """
    m, n = len(ref_seq), len(query_seq)

    M = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Ix = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Iy = np.full((m + 1, n + 1), -np.inf, dtype=float)

    if semi_global:
        M[0, 0] = 0.0
        Ix[0, 0] = -np.inf
        Iy[0, 0] = -np.inf
        for i in range(1, m + 1):
            M[i, 0] = 0.0
            Ix[i, 0] = -np.inf
            Iy[i, 0] = 0.0
        for j in range(1, n + 1):
            M[0, j] = 0.0
            Ix[0, j] = 0.0
            Iy[0, j] = -np.inf
    if semi_global == 'anchor5':
        # 锚定5'端模式: query必须从ref[0]开始，ref不能有前置gap
        # 但允许ref有后置gap（短序列末端gap免费）
        # 用于CARLIN内部区域比对：确保prefix锚定在正确位置
        M[0, 0] = 0.0
        Ix[0, 0] = -np.inf
        Iy[0, 0] = -np.inf
        for i in range(1, m + 1):
            M[i, 0] = -np.inf
            Ix[i, 0] = -np.inf
            Iy[i, 0] = gap_open_profile[i-1] + gap_extend_profile[i-1] * (i-1)
        for j in range(1, n + 1):
            M[0, j] = -np.inf
            Ix[0, j] = -np.inf
            Iy[0, j] = -np.inf
    elif semi_global:
        M[0, 0] = 0.0
        Ix[0, 0] = gap_open_profile[0]
        Iy[0, 0] = gap_open_profile[0]
        for i in range(1, m + 1):
            M[i, 0] = -np.inf
            Ix[i, 0] = -np.inf
            Iy[i, 0] = gap_open_profile[0] + gap_extend_profile[0] * (i - 1)
        for j in range(1, n + 1):
            M[0, j] = -np.inf
            Ix[0, j] = gap_open_profile[0] + gap_extend_profile[0] * (j - 1)
            Iy[0, j] = -np.inf

    # 填充DP矩阵 — 使用位置依赖的gap惩罚 + 行预取加速
    for i in range(1, m + 1):
        Mi_1 = M[i - 1]; Mi = M[i]
        Ixi_1 = Ix[i - 1]; Ixi = Ix[i]
        Iyi_1 = Iy[i - 1]; Iyi = Iy[i]
        r_char = ref_seq[i - 1]

        gp_open_i_1 = gap_open_profile[i - 1]
        gp_extend_i_1 = gap_extend_profile[i - 1]
        gp_open_i = gap_open_profile[min(i, m - 1)]
        gp_extend_i = gap_extend_profile[min(i, m - 1)]
        go_ge_i = gp_open_i + gp_extend_i
        go_ge_i_1 = gp_open_i_1 + gp_extend_i_1

        for j in range(1, n + 1):
            s = match_score if r_char == query_seq[j - 1] else mismatch_penalty
            # M
            a, b, c = Mi_1[j - 1], Ixi_1[j - 1], Iyi_1[j - 1]
            Mi[j] = s + (a if a >= b and a >= c else (b if b >= c else c))
            # Ix
            a, b, c = Mi[j - 1] + go_ge_i, Ixi[j - 1] + gp_extend_i, Iyi[j - 1] + go_ge_i
            Ixi[j] = a if a >= b and a >= c else (b if b >= c else c)
            # Iy
            a, b, c = Mi_1[j] + go_ge_i_1, Ixi_1[j] + go_ge_i_1, Iyi_1[j] + gp_extend_i_1
            Iyi[j] = a if a >= b and a >= c else (b if b >= c else c)

    # 确定最终得分和回溯起点
    if semi_global == 'anchor5':
        # 锚定5'端: query必须耗尽，但ref允许后置gap
        max_score = -np.inf
        max_i, max_j, max_state = m, n, 'M'
        for state, val in [('M', M[m, n]), ('Ix', Ix[m, n]), ('Iy', Iy[m, n])]:
            if val > max_score:
                max_score = val
                max_i, max_j, max_state = m, n, state
    elif semi_global:
        max_score = -np.inf
        max_i, max_j, max_state = 0, 0, 'M'
        for i in range(m + 1):
            for state, val in [('M', M[i, n]), ('Ix', Ix[i, n]), ('Iy', Iy[i, n])]:
                if val > max_score:
                    max_score = val
                    max_i, max_j, max_state = i, n, state
        for j in range(n + 1):
            for state, val in [('M', M[m, j]), ('Ix', Ix[m, j]), ('Iy', Iy[m, j])]:
                if val > max_score:
                    max_score = val
                    max_i, max_j, max_state = m, j, state
    else:
        max_score = max(M[m, n], Ix[m, n], Iy[m, n])
        max_i, max_j = m, n
        max_state = 'M' if M[m, n] >= max(Ix[m, n], Iy[m, n]) else \
                    ('Ix' if Ix[m, n] >= Iy[m, n] else 'Iy')

    # 回溯 — 需要重建 gap_open/extend 以匹配forward pass
    i, j = max_i, max_j
    state = max_state
    aligned_ref, aligned_query = [], []

    while (i > 0 or j > 0) and (semi_global or (i > 0 and j > 0)):
        if state == 'M':
            if i > 0 and j > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append(query_seq[j - 1])
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
                j -= 1
            elif i > 0:
                # ref延伸超出query → 作为deletion
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                i -= 1
                state = 'Iy'
            else:
                break

        elif state == 'Ix':
            if j > 0:
                aligned_ref.append('-')
                aligned_query.append(query_seq[j - 1])
                ri = min(i, m - 1)
                scores = [
                    M[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Ix[i, j - 1] + gap_extend_profile[ri],
                    Iy[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri]
                ]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                j -= 1
            else:
                break

        else:  # Iy
            if i > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                ri = i - 1
                scores = [
                    M[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Ix[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Iy[i - 1, j] + gap_extend_profile[ri]
                ]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
            else:
                break

    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))

    stats = calculate_alignment_stats(aligned_ref, aligned_query)
    stats['score'] = max_score

    return max_score, aligned_ref, aligned_query, stats


def convert_dense_mismatch_to_indel(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    query_seq: str,
    threshold: float = 0.34,
    min_region: int = 3
) -> Tuple[str, str, bool]:
    """
    将高密度mismatch区域转换为insertion (ref侧开gap)。

    当连续比对区域中 mismatch 比例超过 threshold 时，
    认为该query片段实际为insertion事件，将ref对应碱基替换为gap。

    使用较大的滑动窗口和最小mismatch计数来避免误判孤立点突变。

    返回: (corrected_ref, corrected_query, was_modified)
    """
    alen = len(aligned_ref)
    if alen == 0:
        return aligned_ref, aligned_query, False

    # 标记每个位置是否为mismatch (两个都是碱基但不相等)
    is_mismatch = [
        (ar != '-' and aq != '-' and ar != aq)
        for ar, aq in zip(aligned_ref, aligned_query)
    ]

    # 使用较大窗口检测高密度区域 (默认6bp)
    window_size = max(6, min_region)
    half = window_size // 2
    in_dense = False
    dense_regions = []
    region_start = 0

    for i in range(alen):
        start = max(0, i - half)
        end = min(alen, i + half + 1)
        n_bases = sum(1 for k in range(start, end)
                      if aligned_ref[k] != '-' and aligned_query[k] != '-')
        n_mismatch = sum(is_mismatch[start:end])
        # 需要同时满足:
        #   1) 密度 > threshold
        #   2) 至少有2个mismatch (避免孤立点突变被转换)
        density = n_mismatch / n_bases if n_bases > 0 else 0
        is_dense = density > threshold and n_mismatch >= 2

        if is_dense:
            if not in_dense:
                in_dense = True
                region_start = i
        else:
            if in_dense:
                # 只有区域长度 >= min_region 才转换
                if i - region_start >= min_region:
                    dense_regions.append((region_start, i - 1))
                in_dense = False

    if in_dense:
        if alen - region_start >= min_region:
            dense_regions.append((region_start, alen - 1))

    if not dense_regions:
        return aligned_ref, aligned_query, False

    # 执行转换: 在dense region中将ref碱基替换为gap
    ref_list = list(aligned_ref)
    qry_list = list(aligned_query)
    was_modified = False

    for start, end in dense_regions:
        for k in range(start, end + 1):
            if ref_list[k] != '-':
                ref_list[k] = '-'
                was_modified = True

    return ''.join(ref_list), ''.join(qry_list), was_modified


def filter_point_mutations(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    cutsites: List[CutsiteRegion],
    window: int = 3
) -> Tuple[str, str, int]:
    """
    过滤点突变: cutsites ±window 范围外的错配矫正为ref。

    规则:
      1. cutsite ±window 内的点突变 → 保留 (可能是真实编辑)
      2. 紧邻gap的点突变 → 保留 (indel边界处的可能是真实突变)
      3. 其余点突变 → 矫正为ref碱基 (视为假阳性/测序错误)

    返回: (corrected_ref, corrected_query, n_corrected)
    """
    if len(aligned_ref) != len(aligned_query):
        return aligned_ref, aligned_query, 0

    result_query = list(aligned_query)
    n_corrected = 0
    ref_pos = 0  # 参考序列上的坐标 (不含gap)

    for i, (ar, aq) in enumerate(zip(aligned_ref, aligned_query)):
        if ar == '-' or aq == '-':
            if ar != '-':
                ref_pos += 1
            continue

        if ar != aq:
            # 检查是否在 cutsite ±window 内
            in_window = any(cs.start - window <= ref_pos <= cs.end + window for cs in cutsites)
            if not in_window:
                # 检查是否紧邻gap
                adj_to_gap = False
                if i > 0 and (aligned_ref[i - 1] == '-' or aligned_query[i - 1] == '-'):
                    adj_to_gap = True
                if i < len(aligned_ref) - 1 and (aligned_ref[i + 1] == '-' or aligned_query[i + 1] == '-'):
                    adj_to_gap = True

                if not adj_to_gap:
                    result_query[i] = ar
                    n_corrected += 1

        if ar != '-':
            ref_pos += 1

    return aligned_ref, ''.join(result_query), n_corrected


def correct_repetitive_misalignment(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
) -> Tuple[str, str, bool]:
    """
    矫正因重复序列导致的错误比对。

    当参考序列中同一段序列出现在多个位置时，DP算法可能将query片段
    匹配到远端的错误副本，而让近端正确的副本位置留下空缺（deletion）。
    此函数检测该模式并将query碱基从错误位置搬迁到正确位置。

    当前处理的重复序列:
      - ACTGCACGACAGTCG (Target1 ↔ Target9, 15bp)
      - ACTCGCG (Target2 ↔ Target7, 7bp)
      - ACAGTCG (Target1 ↔ Target3 ↔ Target9, 7bp)
      - GAGCGC (Target4 ↔ Target6, 6bp)
      - GCGACT (Target4 ↔ Target7, 6bp)
      - GATACG (Target5 ↔ Target10, 6bp)
      - ACGCAC (Target7 ↔ Target10, 6bp)
      - CGCGCA (Target2 ↔ Target5, 6bp)
      - CGACTA (Target4 ↔ Target9, 6bp)
      - GACGA (Target1 ↔ Target3, 5bp)

    工作原理:
      1. 扫描aligned_ref找到重复序列的所有副本位置
      2. 检测远端副本是否有query匹配，近端副本是否有空缺
      3. 验证近端target区域并非真实删除（≥40%碱基存在）
      4. 将query碱基从远端搬迁到近端
      5. 处理片段开头的孤立碱基（如 GACTGCACGACAGTCG 中的 G）

    返回: (corrected_ref, corrected_query, was_modified)
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False

    ar, aq = aligned_ref, aligned_query
    alen = len(ar)

    # 构建参考序列位置映射
    pos_map = []
    cur = 0
    for c in ar:
        pos_map.append(cur if c != '-' else -1)
        if c != '-':
            cur += 1

    # Target区域定义（用于检查正确位置附近是否有query碱基）
    target_regions = [
        (28, 48), (55, 75), (82, 102), (109, 129), (136, 156),
        (163, 183), (190, 210), (217, 237), (244, 264), (271, 291)
    ]

    def _adjacent_has_bases(ref_pos: int, seg_len: int) -> bool:
        """检查正确位置±3范围内是否有query碱基，避免搬入完全删除的区域"""
        for offset in range(ref_pos - 3, ref_pos + seg_len + 3):
            if 0 <= offset < len(ref_seq):
                col = next((ci for ci, p in enumerate(pos_map) if p == offset), None)
                if col is not None and aq[col] != '-':
                    return True
        return False

    # (repeat_sequence, specific_pairs) 配置
    #   specific_pairs=None  → 自动查找所有副本，第一个为正确位置
    #   specific_pairs=[...] → 使用指定的(正确位置, 错误位置)对
    repeats_config = [
        ("ACTGCACGACAGTCG", None),   # T1-T9
        ("ACTCGCG", None),            # T2-T7
        ("GAGCGC", None),             # T4-T6
        ("GCGACT", None),             # T4-T7
        ("GATACG", None),             # T5-T10
        ("ACGCAC", None),             # T7-T10
        ("CGCGCA", None),             # T2-T5
        ("CGACTA", [(109, 258)]),      # T4-T9 (ref[81]=CGACTA is normal T3上下文，不应作为correct位置)
        ("ACAGTCG", [(37, 86), (37, 253), (86, 253)]),  # T1-T3, T1-T9, T3-T9
        ("GACGA", None),               # T1-T3 (5bp at ref[43]和ref[97])
        ("ACTA", [(260, 83), (260, 125), (125, 83)]),  # T9-T3, T9-T4, T4-T3 (4bp)
    ]

    new_aq = list(aq)
    overall_modified = False

    for repeat, specific_pairs in repeats_config:
        rlen = len(repeat)
        if rlen > alen:
            continue

        # 生成 (correct_rp, wrong_rp) 对
        pairs = []
        if specific_pairs is not None:
            pairs.extend(specific_pairs)
        else:
            copies = []
            p = 0
            while True:
                p = ref_seq.find(repeat, p)
                if p < 0:
                    break
                copies.append(p)
                p += 1
            if len(copies) < 2:
                continue
            for wp in copies[1:]:
                pairs.append((copies[0], wp))

        for correct_rp, wrong_rp in pairs:
            # 通过pos_map逐碱基定位，解决ar列中有gap的情况（query插入导致偏移）
            wrong_cols = [next((ci for ci, p in enumerate(pos_map)
                                if p == wrong_rp + k), -1) for k in range(rlen)]
            correct_cols = [next((ci for ci, p in enumerate(pos_map)
                                  if p == correct_rp + k), -1) for k in range(rlen)]

            if any(c < 0 for c in wrong_cols) or any(c < 0 for c in correct_cols):
                continue

            # 验证ar在错误副本位置确实是该重复序列
            if not all(ar[wrong_cols[k]] == repeat[k] for k in range(rlen)):
                continue

            # 检查query在错误位置是否有匹配（使用原始aq）
            match_cnt = sum(1 for k in range(rlen)
                            if aq[wrong_cols[k]] == repeat[k]
                            and aq[wrong_cols[k]] != '-')
            if match_cnt < rlen * 0.5:
                continue

            # 检查正确位置是否仍有空缺（使用new_aq，之前的搬迁可能已填充）
            gaps_at_correct = sum(1 for k in range(rlen)
                                  if new_aq[correct_cols[k]] == '-')
            if gaps_at_correct < rlen * 0.5:
                continue

            # 验证正确位置前后有query碱基，避免搬入完全删除的区域
            if not _adjacent_has_bases(correct_rp, rlen):
                continue

            # ── 执行搬迁 ──
            modified = False
            for k in range(rlen):
                if new_aq[wrong_cols[k]] != '-' and new_aq[correct_cols[k]] == '-':
                    new_aq[correct_cols[k]] = new_aq[wrong_cols[k]]
                    new_aq[wrong_cols[k]] = '-'
                    modified = True

            # 处理片段开头的孤立G（仅ACTGCACGACAGTCG的前导G需要此步骤）
            if modified and repeat == "ACTGCACGACAGTCG":
                wrong_start_col = wrong_cols[0]
                g_src = None
                for sc in range(wrong_start_col - 1, max(0, wrong_start_col - 30) - 1, -1):
                    if aq[sc] != '-' and aq[sc] == 'G':
                        if ar[sc] == 'G' or ar[sc] == '-':
                            if pos_map[sc] < 0 or abs(pos_map[sc] - correct_rp) >= rlen:
                                if g_src is None:
                                    g_src = sc
                if g_src is not None and correct_rp > 0 and ref_seq[correct_rp - 1] == 'G':
                    if new_aq[correct_cols[0] - 1] == '-':
                        new_aq[correct_cols[0] - 1] = aq[g_src]
                        new_aq[g_src] = '-'
                        modified = True

            if modified:
                overall_modified = True

    if overall_modified:
        return ar, ''.join(new_aq), True
    return aligned_ref, aligned_query, False


def correct_target_misalignments(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
) -> Tuple[str, str, bool]:
    """
    矫正较小片段的跨Target错误比对。

    两类情况:
      1. TAGTAT (ref[219:225], Target8内) 在Target9区域出现时迁回Target8
      2. 单个A碱基在Target9倒数第4位(ref[260])以insertion形式出现时迁回
         Target1倒数第4位(ref[44])

    这些片段出现在错误Target是因为它们与正确Target的序列完全相同，
    DP算法会随机选择匹配位置。
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False

    ar, aq = aligned_ref, aligned_query
    alen = len(ar)

    pos_map = []; cur = 0
    for c in ar:
        pos_map.append(cur if c != '-' else -1)
        if c != '-':
            cur += 1

    new_aq = list(aq)
    modified = False

    # ── 1. TAGTAT → Target8 ──
    tag_seq = 'TAGTAT'
    t_len = 6
    t8_tag_start = ref_seq.find(tag_seq)  # ref[219]
    if t8_tag_start >= 0:
        t8_cs = next((ci for ci in range(alen) if pos_map[ci] == t8_tag_start), None)
        t8_ce = next((ci for ci in range(alen) if pos_map[ci] == t8_tag_start + t_len - 1), None)
        if t8_cs is not None and t8_ce is not None:
            gaps = sum(1 for k in range(t8_cs, t8_ce + 1) if new_aq[k] == '-')
            if gaps >= t_len - 2:
                t8_end_col = next((ci for ci in range(alen) if pos_map[ci] == 236), None)
                if t8_end_col is not None:
                    srch_start = t8_end_col + 1
                    srch_end = min(alen, srch_start + 50)
                    qb = [(ci, new_aq[ci]) for ci in range(srch_start, srch_end)
                          if new_aq[ci] != '-']
                    cols = None
                    for wi in range(len(qb) - t_len + 1):
                        mc = sum(1 for k in range(t_len) if qb[wi + k][1] == tag_seq[k])
                        if mc >= t_len - 1:
                            cols = [qb[wi + k][0] for k in range(t_len)]
                            break
                    if cols:
                        for k in range(t_len):
                            dst = t8_cs + k
                            src = cols[k]
                            if new_aq[src] == tag_seq[k]:
                                new_aq[dst] = new_aq[src]
                                new_aq[src] = '-'
                                modified = True

    # ── 2. 单个A: Target9(ref[260]) → Target1(ref[44]) ──
    a_src_col = next((ci for ci in range(alen) if pos_map[ci] == 260), None)
    a_dst_col = next((ci for ci in range(alen) if pos_map[ci] == 44), None)
    if a_dst_col is not None and a_src_col is not None and new_aq[a_dst_col] == '-':
        srch_start = max(0, a_src_col - 3)
        srch_end = min(alen, a_src_col + 4)
        for ci in range(srch_start, srch_end):
            if ar[ci] == '-' and new_aq[ci] == 'A':
                new_aq[a_dst_col] = 'A'
                new_aq[ci] = '-'
                modified = True
                break
            elif ar[ci] == 'A' and new_aq[ci] == 'A':
                lg = ci > 0 and (ar[ci-1] == '-' or new_aq[ci-1] == '-')
                rg = ci < alen - 1 and (ar[ci+1] == '-' or new_aq[ci+1] == '-')
                if lg and rg:
                    new_aq[a_dst_col] = 'A'
                    new_aq[ci] = '-'
                    modified = True
                    break

    if modified:
        return ar, ''.join(new_aq), True
    return aligned_ref, aligned_query, False


def remove_isolated_matches(
    aligned_ref: str,
    aligned_query: str,
) -> Tuple[str, str, bool]:
    """
    将比对中孤立匹配的单碱基转为gap（deletion）。

    当一个大片段删除区域中间出现一个孤立的单碱基匹配时，
    该匹配会将一个连续deletion切成两个片段。将此孤立匹配转为gap后，
    两侧deletion片段合并为连续删除，使比对结果更简洁且更符合生物学实际。

    孤立匹配判定条件:
      - ref和query在该列有相同碱基（match）
      - 左右相邻列中至少有一列存在gap（ref侧或query侧）
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False

    ar, aq = aligned_ref, aligned_query
    alen = len(ar)
    new_aq = list(aq)
    modified = False

    for i in range(1, alen - 1):
        if ar[i] != '-' and new_aq[i] != '-' and ar[i] == new_aq[i]:
            left_gap = ar[i - 1] == '-' or new_aq[i - 1] == '-'
            right_gap = ar[i + 1] == '-' or new_aq[i + 1] == '-'
            if left_gap and right_gap:
                new_aq[i] = '-'
                modified = True

    if modified:
        return ar, ''.join(new_aq), modified
    return aligned_ref, aligned_query, False


def lineage_tracer_align(
    ref_seq: str,
    query_seq: str,
    cutsites: List[CutsiteRegion],
    match_score: float = 2.0,
    mismatch_penalty: float = -3.0,
    base_gap_open: float = -2.0,
    base_gap_extend: float = -0.1,
    semi_global: bool = True,
    cutsite_gap_scale: float = 1.0,
    flank_gap_scale: float = 2.0,
    far_gap_scale: float = 6.0,
    flank_width: int = 3,
    mismatch_density_threshold: float = 0.34,
    mutation_window: int = 3,
    no_gap_prefix: int = 0
) -> Tuple[float, str, str, Dict]:
    """
    谱系示踪实验专用比对管线。

    整合了结构感知的位置依赖gap惩罚、高密度mismatch→indel转换、
    以及区域感知的点突变过滤, 确保突变的真实性和比对的准确性。

    处理流程:
      1. 构建位置依赖的gap惩罚数组 (cutsite处最低)
      2. 执行位置感知的DP比对
      3. 高密度mismatch区域 → 转换为insertion
      4. 区域感知的点突变过滤
      5. 重新计分统计

    参数:
        ref_seq: 参考序列
        query_seq: 查询序列
        cutsites: CutsiteRegion列表 (需与参考序列坐标对应)
        match_score: 匹配得分 (默认2.0)
        mismatch_penalty: 错配惩罚 (默认-3.0)
        base_gap_open: 基准gap开启惩罚 (默认-2.0)
        base_gap_extend: 基准gap延伸惩罚 (默认-0.1)
        semi_global: 是否半全局比对
        cutsite_gap_scale: cutsite内gap惩罚倍率 (越小gap越容易)
        flank_gap_scale: cutsite侧翼gap惩罚倍率
        far_gap_scale: 远离cutsite区域gap惩罚倍率
        flank_width: cutsite侧翼范围(bp)
        mismatch_density_threshold: mismatch密度阈值 (默认33%)
        mutation_window: 保留点突变的窗口半径 (默认±3bp)
        no_gap_prefix: 前缀长度，此区域内禁止gap开启 (默认0)

    返回: (score, aligned_ref, aligned_query, stats)
    """
    # Step 1: 构建gap惩罚数组
    gap_open_profile, gap_extend_profile = build_gap_penalty_profile(
        ref_length=len(ref_seq),
        cutsites=cutsites,
        base_gap_open=base_gap_open,
        base_gap_extend=base_gap_extend,
        cutsite_scale=cutsite_gap_scale,
        flank_scale=flank_gap_scale,
        far_scale=far_gap_scale,
        flank_width=flank_width
    )

    # Step 1b: 强制前缀区域禁止gap（用于anchor5模式，防止通过级联gap跳过ref起始区域）
    if no_gap_prefix > 0:
        gap_open_profile[:no_gap_prefix] = -1e6
        gap_extend_profile[:no_gap_prefix] = -1e6

    # Step 2: 位置感知比对
    score, aligned_ref, aligned_query, stats = affine_gap_alignment_position_aware(
        ref_seq, query_seq,
        gap_open_profile, gap_extend_profile,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        semi_global=semi_global
    )

    # Step 3: 高密度mismatch → indel 转换
    conv_ref, conv_query, was_modified = convert_dense_mismatch_to_indel(
        aligned_ref, aligned_query, ref_seq, query_seq,
        threshold=mismatch_density_threshold
    )

    # Step 4: 区域感知点突变过滤
    filtered_ref, filtered_query, n_corrected = filter_point_mutations(
        conv_ref if was_modified else aligned_ref,
        conv_query if was_modified else aligned_query,
        ref_seq, cutsites,
        window=mutation_window
    )

    # Step 5: 将孤立匹配转为gap（合并被切断的deletion片段）
    # 在 _align_single 中 correct_repetitive_misalignment 之后执行，
    # 避免在矫正之前移除矫正所需的碱基（如ref[28]的G）
    final_ref, final_query = filtered_ref, filtered_query
    stats['isolated_bases_removed'] = 0
    final_stats = calculate_alignment_stats(final_ref, final_query)
    final_stats['score'] = score
    final_stats['n_mutations_corrected'] = n_corrected
    final_stats['dense_regions_converted'] = was_modified
    final_stats['isolated_bases_consolidated'] = False

    return score, final_ref, final_query, final_stats


def get_amplicon_structure(ref_seq: str) -> List[CutsiteRegion]:
    """
    根据CARLIN扩增子结构自动推断cutsite位置。

    适用于 332bp 标准CARLIN结构:
      Primer5(23) + prefix(5) + 10×[Target(13+7)] + 9×PAM_Linker(7) + postfix(8) + Primer3(33)

    也支持非标准长度参考序列（如去掉Primer的截短版本），
    通过检测 GAGTCG 核心motif的27bp周期来自动推断位置。

    返回: CutsiteRegion列表 (0-indexed, inclusive)
    """
    # 标准结构参数
    target_size = 20        # 每个Target长度
    linker_size = 7         # PAM_Linker长度
    period = target_size + linker_size  # 27bp
    cutsite_offset = 13     # Target内cutsite起始偏移
    cutsite_len = 7         # cutsite长度
    n_targets = 10

    expected_len = 23 + 5 + n_targets * target_size + (n_targets - 1) * linker_size + 8 + 33

    if len(ref_seq) == expected_len:
        # 标准332bp结构 — 直接用公式计算
        target_start = 23 + 5
        cutsites = []
        for i in range(n_targets):
            t_start = target_start + i * period
            cs_start = t_start + cutsite_offset
            cs_end = cs_start + cutsite_len - 1
            cutsites.append(CutsiteRegion(
                name=f"Target{i+1}",
                start=cs_start, end=cs_end
            ))
        return cutsites

    # 非标准长度：通过 GAGTCG motif 周期检测自动推断
    import re
    motif = 'GAGTCG'
    positions = [m.start() for m in re.finditer(motif, ref_seq)]

    if len(positions) < 3:
        print(f"  警告: 无法在参考序列中检测到足够的GAGTCG motif（仅发现{len(positions)}处）")
        print(f"  请手动提供cutsite位置")
        return []

    # 计算相邻motif之间的间隔，推断周期
    diffs = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
    median_diff = sorted(diffs)[len(diffs)//2]  # 中位数间隔

    # 从第一个检测到的motif位置开始，以周期向前/后推断所有cutsite
    first_pos = positions[0]
    # 向前补齐到最近的周期起点
    start = first_pos % period
    while start + period <= first_pos:
        start += period
    while start > first_pos:
        start -= period
    if start < 0:
        start = 0

    cutsites = []
    for i in range(n_targets):
        cs_start = start + i * period
        cs_end = cs_start + cutsite_len - 1
        if cs_start >= len(ref_seq):
            break  # 超出参考序列范围
        if cs_end >= len(ref_seq):
            cs_end = len(ref_seq) - 1
        cutsites.append(CutsiteRegion(
            name=f"Target{i+1}",
            start=cs_start, end=cs_end
        ))

    if cutsites:
        print(f"  auto检测: 发现 {len(cutsites)} 个cutsite区域 (周期={period}bp, motif起始={start})")

    return cutsites


def print_lineage_tracer_alignment(score, aligned_ref, aligned_query, stats,
                                    cutsites, ref_seq, width=100):
    """打印谱系示踪比对结果，标注cutsite区域"""
    print("\n谱系示踪比对结果:")
    print("Cutsite区域标注: [=====]")

    ref_pos = 0
    pos_to_cutsite = {}
    for cs in cutsites:
        for p in range(cs.start, cs.end + 1):
            pos_to_cutsite[p] = cs.name

    for block_start in range(0, len(aligned_ref), width):
        block_end = min(block_start + width, len(aligned_ref))
        ref_slice = aligned_ref[block_start:block_end]
        qry_slice = aligned_query[block_start:block_end]

        # 构建标注行
        annotation = ""
        for k, (r, q) in enumerate(zip(ref_slice, qry_slice)):
            if r == '-' or q == '-':
                annotation += ' '
            elif r == q:
                annotation += '|'
            else:
                annotation += '.'

        print(f"Ref  {block_start:4d}: {ref_slice}")
        print(f"            {annotation}")
        print(f"Query{block_start:4d}: {qry_slice}")
        print()


def print_lineage_stats(stats: Dict, cutsites: List[CutsiteRegion]):
    """打印谱系示踪统计信息"""
    print(f"\n谱系示踪比对统计:")
    print(f"  比对得分: {stats.get('score', 0):.2f}")
    print(f"  比对长度: {stats.get('alignment_length', 0)}")
    print(f"  匹配数: {stats.get('matches', 0)}")
    print(f"  错配数: {stats.get('mismatches', 0)}")
    print(f"  相似度: {stats.get('similarity', 0):.2%}")
    print(f"  一致性: {stats.get('identity', 0):.2%}")
    print(f"  Query中gap数(ref deletion): {stats.get('gaps_in_query', 0)}")
    print(f"  Ref中gap数(query insertion): {stats.get('gaps_in_ref', 0)}")
    print(f"  矫正点突变数: {stats.get('n_mutations_corrected', 0)}")
    print(f"  高密度区域转换: {'是' if stats.get('dense_regions_converted', False) else '否'}")

    if cutsites:
        print(f"\n  Cutsite区域 ({len(cutsites)}个):")
        for cs in cutsites:
            print(f"    {cs.name}: {cs.start}-{cs.end}")


def main():
    """
    主函数：示例使用 — 展示标准比对和谱系示踪比对
    """
    # 示例序列（用户提供的CARLIN序列）
    reference = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )

    # 带大片段缺失的模拟query: 删除了Target3-7区域
    query_del = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )
    # 在cutsite处有突变+小插入的模拟query
    query_mut = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATAGCTACGTGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )

    print("=" * 80)
    print("谱系示踪比对算法示例")
    print("=" * 80)

    # ===== 获取结构信息 =====
    cutsites = get_amplicon_structure(reference)
    if cutsites:
        print(f"\n检测到 {len(cutsites)} 个cutsite区域:")
        print(f"  坐标区间: {cutsites[0].start}-{cutsites[-1].end}")
        print(f"  每个cutsite: 7bp (位置 {cutsites[0].start} 到 {cutsites[-1].end})")

    # ===== 示例1: 谱系示踪比对（带大片段缺失）=====
    print("\n" + "-" * 80)
    print("示例1: 谱系示踪比对 - 大片段缺失 (Target3-7)")
    print("-" * 80)

    print(f"\n序列信息:")
    print(f"  Reference长度: {len(reference)} bp")
    print(f"  Query长度: {len(query_del)} bp")
    print(f"  缺失量: {len(reference) - len(query_del)} bp")

    score1, ar1, aq1, stats1 = lineage_tracer_align(
        reference, query_del, cutsites,
        match_score=2.0, mismatch_penalty=-3.0,
        base_gap_open=-2.0, base_gap_extend=-0.1,
        cutsite_gap_scale=1.0,
        flank_gap_scale=2.0,
        far_gap_scale=6.0,
        flank_width=3
    )

    print_lineage_stats(stats1, cutsites)
    print_alignment(ar1, aq1, width=80)

    # ===== 示例2: 谱系示踪比对（带点突变 + 插入）=====
    print("\n" + "-" * 80)
    print("示例2: 谱系示踪比对 - cutsite处突变 (+ 保守区假阳性突变)")
    print("-" * 80)

    print(f"\n序列信息:")
    print(f"  Reference长度: {len(reference)} bp")
    print(f"  Query长度: {len(query_mut)} bp")

    score2, ar2, aq2, stats2 = lineage_tracer_align(
        reference, query_mut, cutsites,
        match_score=2.0, mismatch_penalty=-3.0,
        base_gap_open=-2.0, base_gap_extend=-0.1,
        cutsite_gap_scale=1.0,
        flank_gap_scale=2.0,
        far_gap_scale=6.0,
        flank_width=3,
        mismatch_density_threshold=0.33,
        mutation_window=3
    )

    print_lineage_stats(stats2, cutsites)
    print_alignment(ar2, aq2, width=80)

    # ===== 对比: 标准比对 vs 谱系示踪比对 =====
    print("\n" + "=" * 80)
    print("对比分析: 标准算法 vs 谱系示踪算法")
    print("=" * 80)

    for q_name, query_seq in [("缺失型", query_del), ("突变型", query_mut)]:
        print(f"\n--- {q_name}序列 ---")

        # 标准算法
        s_std, ar_std, aq_std, st_std = affine_gap_alignment(
            reference, query_seq,
            match_score=2.0, mismatch_penalty=-3.0,
            gap_open=-2.0, gap_extend=-0.1,
            semi_global=True
        )

        # 谱系示踪算法
        s_lt, ar_lt, aq_lt, st_lt = lineage_tracer_align(
            reference, query_seq, cutsites,
            match_score=2.0, mismatch_penalty=-3.0,
            base_gap_open=-2.0, base_gap_extend=-0.1
        )

        print(f"  标准: 得分={s_std:.2f}, 匹配={st_std['matches']}, 错配={st_std['mismatches']}, "
              f"gap_ref={st_std['gaps_in_ref']}, gap_qry={st_std['gaps_in_query']}, 相似={st_std['similarity']:.2%}")
        print(f"  谱系: 得分={s_lt:.2f}, 匹配={st_lt['matches']}, 错配={st_lt['mismatches']}, "
              f"gap_ref={st_lt['gaps_in_ref']}, gap_qry={st_lt['gaps_in_query']}, 相似={st_lt['similarity']:.2%}")
        print(f"  矫正: {st_lt.get('n_mutations_corrected', 0)} 个假阳性突变被矫正")


if __name__ == "__main__":
    main()