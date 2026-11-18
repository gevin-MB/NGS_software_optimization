# split_barcode — 多线程 FASTQ barcode 拆分工具（C 版）

对 `jinghq/tools` 下 Perl 拆分程序（`MGI_SplitDualBarcodes.pl`、`SplitDualBarcodes_forSE.pl`、
`MGI_SplitBarcode_PE150.pl`）的 C 重写：多线程拆分、支持单端/双端、单/双 barcode、barcode 位置可指定。

## 编译

```bash
make CC=/usr/bin/gcc          # 或 make（conda gcc 缺 -lz 时指定系统 gcc）
```

依赖：zlib、pthread（无其他第三方库）。单文件 `split_barcode.c`，约 900 行。

## 用法

```bash
# 双端双 barcode（MGI_SplitDualBarcodes.pl 完全兼容，布局: col2 在 f, col1 在 f+i+1，均取反向互补）
./split_barcode -r1 R1.fq.gz -r2 R2.fq.gz -f 101 -i 19 -b barcode.list -rc Y -e 1 -t 8 --mgi-dual -o outdir

# 双端双 barcode，显式指定每段位置（1-based），列顺序 = 列表列顺序
./split_barcode -r1 R1.fq.gz -r2 R2.fq.gz -b barcode.list -p "121,101" -rc Y -e 1 -t 8 -o outdir

# 单端单 barcode（SplitDualBarcodes_forSE.pl SE 用法）
./split_barcode -r1 R1.fq.gz -f 5 -b se_barcode.list -rc Y -e 1 -t 8 -o outdir
```

### 参数

| 参数 | 说明 | 默认 |
|---|---|---|
| `-r1/--read1` | read1.fq / .fq.gz（必需） | — |
| `-r2/--read2` | read2.fq / .fq.gz；不提供即为单端 SE | — |
| `-b/--barcodeList` | barcode 列表（必需） | — |
| `-f/--firstCycle` | barcode 起始 cycle（1-based） | — |
| `-i/--interval` | 双 barcode 间隔；传统模式第二段起始 = f+len1+i | — |
| `-p/--positions` | 显式指定各段位置 `"101:10,131:8"` 或 `"101,131"`（len 缺省取列表长度），优先于 -f/-i | — |
| `-e/--errNum` | 容错 mismatch 数（邻居哈希预展开；对每段独立允许 ≤e 个错） | PE=2, SE=1 |
| `-w/--bcRead` | barcode 位于哪条 read | PE=2, SE=1 |
| `-l/--lastCycle` | 单 barcode 结束 cycle，配合 -f 指定区间（兼容 MGI_SplitBarcode.pl） | — |
| `-x/--no-unbarcoded` | 不输出 unbarcoded 文件（仅计数，兼容 MGI_SplitBarcode.pl 行为） | — |
| `--cut-before` | barcode read 输出仅保留 barcode 之前部分（兼容 MGI_SplitBarcode.pl，常与 `-u 9` 同用） | 关 |
| `-u/--umiLen` | barcode 区之后取 N 碱基作为 UMI 附加到 read 名（`#UMI`，保留原 `/x` 后缀） | 0（关闭） |
| `-rc/--revcom` | 列表序列反向互补后匹配 | Y |
| `-c/--compress` | 输出 .fq.gz | Y |
| `-z/--level` | gzip 压缩级别 | 4 |
| `-t/--threads` | 工作线程数（另有每输出文件 1 个压缩写线程） | 4 |
| `-o/--outdir` | 输出目录 | ./ |
| `--mgi-dual` | 双 barcode 布局兼容 MGI_SplitDualBarcodes.pl（col2@f, col1@f+i+1） | 关 |

### barcode.list 格式（空白分隔，# 为注释）

```
#barcodeNum  barcodeSeq1  [barcodeSeq2]
1            ATGCATCTAA   TATAGCCTAG
2            AGCTCTGGAC   CTCTATCGTC
```

首行序列列数决定单/双 barcode 模式；单端时每行 1 列序列。

## 输出

- `前缀_样本_1.fq.gz` / `_2.fq.gz`（SE 为 `前缀_样本.fq.gz`）：barcode 区从序列和质量中切除（两段及其间隔整体切除）
- `前缀_unbarcoded_1/2.fq.gz`：未匹配 reads 原样保留
- `BarcodeStat.txt`：每样本 Correct / Corrected / Total / Pct（格式同 Perl）
- `TagStat.txt`：观测 tag 频数表（按频数降序，格式同 Perl）

## 与 Perl 版的行为差异（已知且有意为之）

1. **容错匹配口径**：C 版对每段 barcode 独立允许 ≤e 个 mismatch（与 usage 文档一致）。
   Perl `MGI_SplitDualBarcodes.pl` 的 `bar_hash` 双 barcode 分支实为"两段各 1 错且替换同一碱基"，
   导致 `-e 1` 时精确匹配 reads 大量漏匹配（见 bench：Perl Correct=0 vs C Correct=1699705）。
   C 版为修正后的语义；使用 `--mgi-dual` 时提取/布局与 Perl 一致，但容错仍为 C 口径。
2. Perl `SplitDualBarcodes_forSE.pl`/`MGI_SplitBarcode_PE150.pl` 的 tag 提取取 `len+1` 字符
   （疑似 off-by-one），C 版取精确长度。
3. 输出 read 的切头/UMI 行为：C 版统一为"切除 barcode 区 + 可选 -u UMI"；MGI Perl 会把 read2
   截断到 barcode 前、双端加 `#UMI/1`，与协议相关，如需完全一致请沿用原 Perl。
4. `-i` 语义：C 传统模式 `第二段起始 = f + len1 + i`；MGI Perl 为 `第二段 1-based 起始 = f + i + 1`。
   兼容 Perl 数据请用 `--mgi-dual`（等价 `-i 19` 对应原 Perl `-i 19`）。
5. **歧义容错邻居**：同一序列落在多个 barcode 的容错范围内时，Perl 后写覆盖先写（取决于列表顺序），
   C 版按未匹配处理（更安全）。
6. TagStat 未匹配标签列：C 版为 `unknown`（MGI_SplitBarcode.pl 拼写为 `unknow`）。

## MGI_SplitBarcode.pl（双端单 barcode）兼容用法

```bash
# 原: perl MGI_SplitBarcode.pl -r1 R1.fq.gz -r2 R2.fq.gz -e 1 -f 101 -l 110 -b barcode.list -o out
./split_barcode -r1 R1.fq.gz -r2 R2.fq.gz -e 1 -f 101 -l 110 -b barcode.list \
    -rc Y -u 9 --cut-before -x -t 8 -o out
```

已在 2 万对合成数据上验证：BarcodeStat 与 Perl 完全一致，8 个样本 fastq 内容逐字节一致
（多线程输出不保序，按内容比对）。

注意 Perl `MGI_SplitBarcode.pl` 两处缺陷，C 版为修正语义：
1. 第 97 行 `$barcode_len=length($tmp[1])` 缺分号，与下一行合并成
   `$barcode_len = length(...) & bar_hash(...)`——UMI 偏移变成 `10 & 引用地址`（0/2/8/10，随内存地址
   变化，非确定性）；正确本意应为 barcode 后取 9 碱基，即 C 版 `-u 9`。
2. read2 截断长度硬编码 100（恰好等于 fc=101 时的 fc-1），fc≠101 时行为不符预期；C 版
   `--cut-before` 按 barcode 实际起始位置截断。

## 性能（200 万对 PE、4 样本、同机对比）

| 程序 | 墙钟 | CPU 时间 | 结果 |
|---|---|---|---|
| `MGI_SplitDualBarcodes.pl -e 1` | 17.1 s | 127 s | 见上 |
| `split_barcode -t 8` | **4.4 s** | **27 s** | BarcodeStat 与 Perl 逐位一致，峰值内存 ~77 MB |

正确性：4 万对 PE（双 barcode，含 1 错配/未知）与 3 万条 SE（单 barcode、rc、含错配）合成数据，
拆分计数与注入真值完全一致（见 `gen_test.py`）。

## 实现

- 读线程 gzread 解压 → 2 万条/块 → 工作队列（上限 8 块在飞）→ N 工作线程匹配（预展开邻居哈希，
  开放寻址）→ 每输出文件独立写线程 zlib 压缩（1MB 块，`gzbuffer` 512KB-1MB）
- 统计合并：每 worker 独立 tag 哈希与 correct/corrected 数组，结束后归并
