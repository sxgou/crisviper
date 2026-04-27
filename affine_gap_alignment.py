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
from typing import Tuple, Dict, List


def affine_gap_alignment(ref_seq: str, query_seq: str, 
                         match_score: float = 2.0, 
                         mismatch_penalty: float = -3.0,
                         gap_open: float = -2.0, 
                         gap_extend: float = -0.1,
                         semi_global: bool = True) -> Tuple[float, str, str, Dict]:
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
    if semi_global:
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
    
    # 填充DP矩阵
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # 匹配/错配得分
            s = match_score if ref_seq[i - 1] == query_seq[j - 1] else mismatch_penalty
            
            # M状态：来自匹配/错配
            M[i, j] = s + max(M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
            
            # Ix状态：gap在ref中（query插入）
            Ix[i, j] = max(
                M[i, j - 1] + gap_open + gap_extend,
                Ix[i, j - 1] + gap_extend,
                Iy[i, j - 1] + gap_open + gap_extend
            )
            
            # Iy状态：gap在query中（ref缺失）
            Iy[i, j] = max(
                M[i - 1, j] + gap_open + gap_extend,
                Ix[i - 1, j] + gap_open + gap_extend,
                Iy[i - 1, j] + gap_extend
            )
    
    # 确定最终得分和回溯起点
    if semi_global:
        # 取最后一行和最后一列的最大值
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
    
    while (i > 0 or j > 0) and (semi_global or (i > 0 and j > 0)):
        if state == 'M':
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append(query_seq[j - 1])
            # 前一个状态
            scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
            prev = np.argmax(scores)
            state = ['M', 'Ix', 'Iy'][prev]
            i -= 1
            j -= 1
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


def main():
    """
    主函数：示例使用
    """
    # 示例序列（用户提供的CARLIN序列）
    reference = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )
    
    query = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACGTATCTACGTGCGTATCTACGTGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )
    
    print("=" * 80)
    print("仿射gap惩罚序列比对算法示例")
    print("=" * 80)
    
    print(f"\n序列信息:")
    print(f"  Reference长度: {len(reference)} bp")
    print(f"  Query长度: {len(query)} bp")
    
    # 运行比对
    print("\n运行比对算法...")
    score, aligned_ref, aligned_query, stats = affine_gap_alignment(
        reference, query,
        match_score=2.0,
        mismatch_penalty=-3.0,
        gap_open=-2.0,
        gap_extend=-0.1,
        semi_global=True
    )
    
    # 打印结果
    print_stats(stats)
    print_alignment(aligned_ref, aligned_query, width=100)
    
    # 示例：不同参数对比
    print("\n" + "=" * 80)
    print("不同参数对比示例")
    print("=" * 80)
    
    param_sets = [
        {"name": "默认参数", "match": 2.0, "mismatch": -3.0, "open": -2.0, "extend": -0.1},
        {"name": "严格匹配", "match": 3.0, "mismatch": -5.0, "open": -2.0, "extend": -0.1},
        {"name": "鼓励插入", "match": 2.0, "mismatch": -3.0, "open": -0.5, "extend": -0.01},
        {"name": "全局比对", "match": 2.0, "mismatch": -3.0, "open": -2.0, "extend": -0.1},
    ]
    
    for params in param_sets:
        print(f"\n{params['name']}:")
        semi_global = params['name'] != "全局比对"
        score, _, _, stats = affine_gap_alignment(
            reference, query,
            match_score=params['match'],
            mismatch_penalty=params['mismatch'],
            gap_open=params['open'],
            gap_extend=params['extend'],
            semi_global=semi_global
        )
        print(f"  得分: {score:.2f}, 匹配: {stats['matches']}, 错配: {stats['mismatches']}, "
              f"相似度: {stats['similarity']:.2%}")


if __name__ == "__main__":
    main()