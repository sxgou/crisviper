# Contributing / 贡献指南

## Reporting Issues / 报告问题

- 搜索已有 Issue 避免重复
- 提供完整的重现步骤、输入数据和预期输出
- 附上 `crisviper --version` 输出

## Pull Requests

1. Fork 仓库并创建功能分支：`git checkout -b feature/my-feature`
2. 确保代码通过 lint：`ruff check crisviper/`
3. 确保全部测试通过：`python -m pytest tests/`
4. 提交 PR 到 `dev` 分支

## Development Setup

```bash
git clone https://github.com/sxgou/crisviper
cd crisviper
python -m venv venv
source venv/bin/activate
pip install -e .
pip install ruff pytest
```

## Code Style

- 遵循 PEP 8
- 所有新函数需有 docstring（类型注解 + 功能说明）
- 所有变更需有对应测试
- 使用 `ruff check crisviper/` 检查代码质量

## Testing

```bash
python -m pytest tests/          # 全部测试
python -m pytest tests/ -v       # 详细输出
python -m pytest tests/ -k test_name  # 运行单个测试
```

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
