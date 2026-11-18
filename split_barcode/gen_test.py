#!/usr/bin/env python3
"""生成模拟 fastq 验证 split_barcode"""
import random, gzip, os

random.seed(42)
OUT = os.path.dirname(os.path.abspath(__file__))
TMP = os.path.join(OUT, "test")
os.makedirs(TMP, exist_ok=True)

def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

# 4 个样本, 双 10bp barcode
barcodes = [
    ("1", "ATGCATCTAA", "TATAGCCTAG"),
    ("2", "AGCTCTGGAC", "CTCTATCGTC"),
    ("3", "TTGACCAGTA", "GGATCAACCT"),
    ("4", "CCATGTGAAG", "TCAAGGCTTA"),
]

F = 101        # barcode1 起点 (1-based)
I = 10         # 间隔
POS1 = F - 1              # 0-based
POS2 = F + 10 + I - 1     # barcode2 起点 0-based

def randseq(n): return "".join(random.choice("ACGT") for _ in range(n))

def make_record(name):
    seq1, qual1 = randseq(150), "".join(random.choice("I") for _ in range(150))
    seq2, qual2 = randseq(150), "".join(random.choice("I") for _ in range(150))
    return seq1, qual1, seq2, qual2

# 统计真值: {sample: (exact, corrected)}, unknown 数
truth = {b[0]: [0, 0] for b in barcodes}
unknown = 0
r1, r2 = [], []
n = 40000
for k in range(n):
    seq1, qual1, seq2, qual2 = make_record(k)
    name = "@read_%d" % k
    r = random.random()
    if r < 0.10:      # 10% 未知 barcode
        bc1, bc2 = randseq(10), randseq(10)
        unknown += 1
    else:
        b = random.choice(barcodes)
        bc1, bc2 = b[1], b[2]
        if r > 0.85:  # 15% 中随机引入 1 个 mismatch(在前段), -e 1 应纠正
            p = random.randrange(10)
            orig = bc1[p]
            bc1 = bc1[:p] + random.choice([c for c in "ACGT" if c != orig]) + bc1[p+1:]
            truth[b[0]][1] += 1
        else:
            truth[b[0]][0] += 1
    seq2 = seq2[:POS1] + bc1 + seq2[POS1+10:POS2] + bc2 + seq2[POS2+10:]
    qual2 = qual2[:POS1] + "I"*20 + qual2[POS1+20:]
    r1.append("%s/1\n%s\n+\n%s\n" % (name, seq1, qual1))
    r2.append("%s/2\n%s\n+\n%s\n" % (name, seq2, qual2))

with open(os.path.join(TMP, "test_1.fq"), "w") as f: f.writelines(r1)
with open(os.path.join(TMP, "test_2.fq"), "w") as f: f.writelines(r2)
with open(os.path.join(TMP, "barcode.list"), "w") as f:
    f.write("#barcodeNum\tbarcodeSeq1\tbarcodeSeq2\n")
    for b in barcodes:
        f.write("%s\t%s\t%s\n" % b)

# SE 数据: 单 10bp barcode 在 read1 第 5 位, rc=Y (read 携带 revcomp)
se_bar = [("S%d" % i, b) for i, b in enumerate(["GGTACCATGG", "TTAGCACTCG", "CCATGATTAC"])]
se_truth = {b[0]: [0, 0] for b in se_bar}
se_unknown = 0
se = []
for k in range(30000):
    seq1, qual1, _, _ = make_record(k)
    name = "@se_%d" % k
    r = random.random()
    if r < 0.10:
        bc = randseq(10); se_unknown += 1
    else:
        b = random.choice(se_bar); bc = b[1]
        if r > 0.85:
            p = random.randrange(10); orig = bc[p]
            bc = bc[:p] + random.choice([c for c in "ACGT" if c != orig]) + bc[p+1:]
            se_truth[b[0]][1] += 1
        else:
            se_truth[b[0]][0] += 1
    bc_r = rc(bc)
    seq1 = seq1[:4] + bc_r + seq1[14:]
    qual1 = qual1[:4] + "I"*10 + qual1[14:]
    se.append("%s\n%s\n+\n%s\n" % (name, seq1, qual1))
with open(os.path.join(TMP, "se.fq"), "w") as f: f.writelines(se)
with open(os.path.join(TMP, "se_barcode.list"), "w") as f:
    for b in se_bar: f.write("%s\t%s\n" % b)

print("PE truth:", truth, "unknown:", unknown, "total:", n)
print("SE truth:", se_truth, "unknown:", se_unknown, "total:", 30000)
