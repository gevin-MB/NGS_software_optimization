# NGS Software Optimization

本仓库用于沉淀 NGS 数据处理流程中自研优化的高性能工具，替代或补充原有 Perl/Shell 实现。每个工具一个独立子目录，互不依赖，均可单独编译部署。

## 工具列表

| 目录 | 工具 | 功能 | 状态 |
|---|---|---|---|
| `split_barcode/` | split_barcode | FastQ barcode 拆分（单/双条码，C 多线程） | 已验证，可替换生产 |

## split_barcode 简介

FastQ 条码拆分工具，多线程 + gzip 写出流水线，兼容华大 MGI 拆分脚本语义：

- 单条码 / 双条码（`--mgi-dual`，兼容 `MGI_SplitDualBarcodes*.pl` 含 NAD 变体）
- `-e` 容错匹配（默认 PE=2/SE=1），N 碱基按一次错配处理
- 精确/容错两级哈希查找，歧义键自动检测并丢弃（不覆盖）

### 编译

```bash
cd split_barcode
make CC=/usr/bin/gcc     # -O3 -march=native -std=gnu99，依赖 zlib / pthread
```

### 性能（对比原 Perl 生产脚本）

| 数据量 | Perl 生产封装 | split_barcode | 提速 |
|---|---|---|---|
| 1.31 亿对 PE150（单 lane） | 8:58 | 5:51（单进程 8 线程） | 1.5× |
| 20.5 亿对 PE100（4 lane 满跑） | ≈90 分钟 | 21:20（4 进程 × 8 线程） | 4.2× |

精确匹配层与 Perl 结果逐样本完全一致（小芯片 1.23 亿条、大芯片 19.8 亿条）；
`-e` 容错层存在 ~0.02% 归属差异，源于 Perl `bar_hash` 的实现缺陷（N 键要求两段同 N、双段变异共用同一替换碱基），C 侧行为更合理，详见各工具目录内的基准报告。

## 新增工具约定

- 每个工具一个子目录：`<tool>/`（源码 + Makefile + README + 测试数据生成脚本）
- 编译产物、fastq、大测试数据一律不入库（见 `.gitignore`）
- 上线前在真实下机数据上与原实现做逐条 read 对账，报告放在工具目录内
