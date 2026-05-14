"""Dataclass definitions for amplicon structure configuration."""

from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class CutsiteRegion:
    """cutsite区域定义 (0-indexed, inclusive)"""
    name: str               # 如 "Target1"
    start: int              # cutsite起始位置 (参考序列上的坐标)
    end: int                # cutsite结束位置 (参考序列上的坐标)


@dataclass
class AmpliconConfig:
    """扩增子结构配置 — 集中管理CARLIN特异性参数。

    所有与特定扩增子设计相关的硬编码参数集中在此，
    便于适配不同的多靶点谱系示踪实验设计。
    """
    primer5_len: int = 23       # 5' 端引物长度 (bp)
    primer3_len: int = 33       # 3' 端引物长度 (bp)
    prefix: str = "CGCCG"       # 5' 前缀序列 (Primer5之后)
    postfix_len: int = 8        # 3' 后缀长度 (Primer3之前)

    target_size: int = 20       # 每个Target长度 (保守区13 + cutsite 7)
    linker_size: int = 7        # PAM_Linker长度
    n_targets: int = 10         # Target数量
    cutsite_offset: int = 13    # Target内cutsite起始偏移
    cutsite_len: int = 7        # cutsite长度

    dual_anchor_tolerance: int = 4  # 双端锚定容许错配数

    @property
    def period(self) -> int:
        """Target + Linker 周期长度"""
        return self.target_size + self.linker_size

    @property
    def expected_full_length(self) -> int:
        """标准扩增子全长 (有Primer)"""
        return (self.primer5_len + len(self.prefix) +
                self.n_targets * self.target_size +
                (self.n_targets - 1) * self.linker_size +
                self.postfix_len + self.primer3_len)

    @classmethod
    def carlin_standard(cls) -> "AmpliconConfig":
        """标准CARLIN扩增子 (332bp) 配置"""
        return cls()

    @classmethod
    def from_json(cls, path: str) -> "AmpliconConfig":
        """从JSON配置文件加载"""
        import json
        with open(path) as f:
            data = json.load(f)
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


