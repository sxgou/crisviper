# COMPLEX → INDEL 重构设计

## 背景

当前代码使用 `COMPLEX` 类型表示相邻插入+删除的复合突变事件。存在以下问题：

1. **INS length bug**：合并计算时，插入长度被当作 ref 碱基跨度，但 INS 不占据 ref 位置
2. **仅两两合并**：只处理相邻的 INS+DEL 对，不支持 INS+DEL+INS 等多事件连续合并
3. **命名歧义**：`COMPLEX` 语义模糊，`INDEL` 更准确地描述事件类型

## 改动概览

- 枚举 `COMPLEX` → `INDEL`
- 重写合并算法为贪心分组，支持任意数量的连续相邻事件
- 将 SUB 也纳入分组合并，只要组内至少有一个 INS 或 DEL
- 删除已废弃的 MATLAB 兼容函数
- 更新所有引用处的类型名

## 合并算法

### 相邻判定

对每个事件定义其在 ref 上的区间（exclusive end）：

| 类型 | start | end (exclusive) |
|------|-------|-----------------|
| DEL  | `ref_pos` | `ref_pos + length` |
| INS  | `ref_pos` | `ref_pos + 1` |
| SUB  | `ref_pos` | `ref_pos + 1` |

两个事件相邻当且仅当它们的闭区间重叠或恰好差 1bp：
```
A.start <= B.end - 1 + 1  即  A.start <= B.end
且 B.start <= A.end
```
等价于 `max(A.start, B.start) <= min(A.end, B.end) + 1`。

### 贪心分组

```
遍历 events（已按 ref_pos 排序）：
  初始 current_group = [events[0]]
  对后续事件 e：
    如果 e 与 current_group 中所有事件相邻（与 group_end 比较即可）：
      current_group.append(e)
    否则：
      消化 current_group → 输出结果
      开始新组 [e]

消化规则：
  组内至少有一个 INS 或 DEL → 合并为 INDEL
  全为 SUB → 逐个输出为独立 SUB
```

### 合并坐标计算

```
ref_start = min(e.ref_pos for e in group)
ref_end   = max(e.end_exclusive for e in group)
length    = ref_end - ref_start

// 示例: DEL(5,3) + INS(6) + SUB(9)
//   ends: 8, 7, 10
//   ref_start=5, ref_end=10, length=5 → covers [5,9]
```

## 波及范围

### 枚举改名

`models.py:21`: `COMPLEX = "complex"` → `INDEL = "indel"`

### mutation.py

| 函数 | 改动 |
|------|------|
| `_merge_adjacent_indels` | 重写为贪心分组算法，修复 INS length bug，纳入 SUB |
| `identify_sequence_events` | 删除（已由 `extract_mutations` 取代） |
| `identify_cas9_events` | 删除（已由 `extract_mutations` + 合并取代） |
| `_find_site_for_bp` | 删除（仅 `identify_cas9_events` 使用） |
| `classify_bp_event` | `COMPLEX` → `INDEL` |
| `annotate_mutation` | `COMPLEX` → `INDEL` |
| `format_mutations_for_display` | `COMPLEX` → `INDEL` |

### caller.py

`_event_structure`: `COMPLEX` → `INDEL`，前缀保持 `"C"`（内部编码，与 INSERTION 的 `"I"` 不冲突）

### reporting.py

`_write_allele_annotations`: `identify_cas9_events` → `extract_mutations`

### __init__.py

`classify_bp_event` 保留并在 `__init__.py` 中继续导出（per-base 分类仍被使用）。
删除 `identify_sequence_events`, `identify_cas9_events`, `_find_site_for_bp` 的导出。

### 测试文件

| 文件 | 改动 |
|------|------|
| `tests/test_mutation.py` | COMPLEX → INDEL；删除 `identify_sequence_events`/`identify_cas9_events` 相关测试 |
| `tests/test_pipeline.py` | COMPLEX → INDEL |
| `tests/test_accuracy.py` | COMPLEX → INDEL |
| `tests/test_gap_exit_bonus.py` | COMPLEX → INDEL |

## 管线顺序（不变）

```
比对 → 背景点突变矫正 → extract_mutations（内部包含 INDEL 合并）→ allele 过滤
```

矫正函数工作在比对序列层面，合并在 MutationEvent 层面，互不影响。

## 测试要点

- `DEL(5,3) + INS(6)` → `INDEL(5, len=3)` 覆盖 [5,7]
- `INS(4) + DEL(5,3)` → `INDEL(4, len=4)` 覆盖 [4,7]
- `DEL(5,3) + INS(8)` → `INDEL(5, len=4)` 覆盖 [5,8]
- `INS(4) + DEL(5,3) + INS(9)` → `INDEL(4, len=6)` 覆盖 [4,9]
- `SUB(2) + INS(4) + DEL(5,3)` → `INDEL(2, len=6)` 覆盖 [2,7]
- `SUB(2) + SUB(5)`（无 INS/DEL）→ 保持独立 SUB 事件
- 间隔 >1bp 的事件不合并
