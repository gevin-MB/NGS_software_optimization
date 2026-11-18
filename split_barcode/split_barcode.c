/*
 * split_barcode.c  —  多线程 FASTQ barcode 拆分工具 (C + pthreads + zlib)
 *
 * Perl 原版: MGI_SplitDualBarcodes.pl / SplitDualBarcodes_forSE.pl
 * 功能:
 *   - 单端(SE) / 双端(PE) FASTQ 拆分, 支持 .fq / .fq.gz 输入
 *   - 单 barcode 或双 barcode (barcode.list 中 1 列或 2 列序列自动识别)
 *   - barcode 位置可指定: -f 起始cycle(1-based) + -i 双barcode间隔,
 *     或 -p 显式指定每段 barcode 的位置(和长度)
 *   - -e 容错 mismatches (邻居哈希预展开), -rc 反向互补, -u UMI 附加到 read 名
 *   - 多线程: 1 读线程 + N 工作线程, gzip 压缩在写线程内完成
 *   - 输出 BarcodeStat.txt / TagStat.txt, 格式与 Perl 版一致
 *
 * 编译: gcc -O3 -march=native -pthread split_barcode.c -o split_barcode -lz
 */
#define _GNU_SOURCE
#include <ctype.h>
#include <errno.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#define MAX_SEG     8
#define MAX_NAME    256
#define MAX_SEQLEN  128
#define OUT_BUF_CAP (1u << 20)

/* ---------------- utils ---------------- */
static void die(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "[split_barcode] ERROR: ");
    va_start(ap, fmt); vfprintf(stderr, fmt, ap); va_end(ap);
    fprintf(stderr, "\n");
    exit(1);
}
static void *xmalloc(size_t n)  { void *p = malloc(n ? n : 1); if (!p) die("out of memory"); return p; }
static void *xcalloc(size_t n, size_t s) { void *p = calloc(n ? n : 1, s); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *p, size_t n) { void *q = realloc(p, n ? n : 1); if (!q) die("out of memory"); return q; }
static char *xstrdup(const char *s) { char *p = strdup(s); if (!p) die("out of memory"); return p; }

static char comp(char c);

static void revcomp(char *s, int len)
{
    for (int i = 0, j = len - 1; i <= j; i++, j--) {
        char a = s[i], b = s[j];
        s[i] = comp(b); s[j] = comp(a);
    }
}
static char comp(char c)
{
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        default:  return 'N';
    }
}

static uint64_t fnv(const char *s, size_t n)
{
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= (unsigned char)s[i]; h *= 1099511628211ULL; }
    return h;
}

/* ---------------- options ---------------- */
static const char *opt_r1 = NULL, *opt_r2 = NULL, *opt_bl = NULL, *opt_outdir = ".";
static int  opt_errnum   = -1;   /* default: PE 2, SE 1 */
static int  opt_first    = 0;    /* -f, 1-based */
static int  opt_interval = -1;   /* -i */
static int  opt_threads  = 4;
static int  opt_level    = 4;
static int  opt_rc       = 1;    /* -rc Y default */
static int  opt_compress = 1;    /* -c Y default */
static int  opt_bcread   = 0;    /* -w 1|2, default: 2 if PE else 1 */
static int  opt_umi      = 0;    /* -u UMI length after barcode span */
static char opt_pos[1024] = "";
static int  opt_mgi = 0;         /* --mgi-dual: 双barcode布局兼容 MGI_SplitDualBarcodes.pl
                                    (col2 在 f, col1 在 f+interval+1) */
static int  opt_last = 0;        /* -l/--lastCycle: 单 barcode 结束 cycle (MGI_SplitBarcode.pl 兼容) */
static int  opt_cut_before = 0;  /* --cut-before: barcode read 输出只保留 barcode 之前部分
                                    (兼容 MGI_SplitBarcode.pl 的 substr(seq2,0,fc-1)) */
static int  opt_no_unbar = 0;    /* -x/--no-unbarcoded: 不输出 unbarcoded 文件(只计数) */
static int  chunk_records = 10000;
static int  max_inflight  = 8;

static const char *g_usage =
"Usage: split_barcode [options]\n"
"  *-r1 --read1 <str>        read1.fq[.gz]\n"
"   -r2 --read2 <str>        read2.fq[.gz]  (不提供则为单端 SE)\n"
"  *-b  --barcodeList <str>  barcode 列表\n"
"   -f  --firstCycle <int>   barcode 起始 cycle (1-based, 默认作用于 read2/SE 的 read1)\n"
"   -i  --interval <int>     双 barcode 时第二段起始 = f + len1 + interval\n"
"   -p  --positions <spec>   显式指定每段位置: \"101:10,131:8\" 或 \"101,131\"\n"
"                            (start 为 1-based; len 缺省取列表序列长度; 给定 -p 时忽略 -f/-i)\n"
"   -e  --errNum <int>       容错 mismatches [默认 PE=2, SE=1]\n"
"   -w  --bcRead <1|2>       barcode 位于哪条 read [默认 PE=2, SE=1]\n"
"   -u  --umiLen <int>       barcode 区之后取 N 碱基作为 UMI 附加到 read 名 (#UMI)\n"
"   -rc --revcom <Y|N>       列表序列取反向互补后匹配 [默认 Y]\n"
"   -c  --compress <Y|N>     输出 gzip 压缩 [默认 Y]\n"
"   -z  --level <int>        gzip 压缩级别 1-9 [默认 4]\n"
"   -l  --lastCycle <int>    单 barcode 结束 cycle (兼容 MGI_SplitBarcode.pl, 需与 -f 同用)\n"
"   -x  --no-unbarcoded      不输出 unbarcoded 文件, 仅计数\n"
"       --cut-before         barcode read 输出仅保留 barcode 之前部分\n"
"                            (兼容 MGI_SplitBarcode.pl 的 read2 截断, 一般与 -u 9 同用)\n"
"   -t  --threads <int>      工作线程数 [默认 4]\n"
"   -o  --outdir <str>       输出目录 [默认 ./]\n"
"   -h  --help\n"
"\n"
"barcode.list 格式 (制表符/空白分隔, # 开头为注释):\n"
"  单 barcode:  #name  seq\n" 
"               1      ATGCATCTAA\n"
"  双 barcode:  #name  seq1  seq2\n"
"               1      ATGCATCTAA  TATAGCCTAG\n";

static int match_opt(const char *arg, const char *sh, const char *lo)
{
    return (sh && strcmp(arg, sh) == 0) || (lo && strcmp(arg, lo) == 0);
}

/* ---------------- barcode data ---------------- */
typedef struct {
    char name[MAX_NAME];
    char seg[MAX_SEG][MAX_SEQLEN];
    int  slen[MAX_SEG];
} Barcode;

static Barcode *g_bc   = NULL;
static int      g_nbc  = 0;
static int      g_nseg = 0;

typedef struct { int start; int len; } SegPos;   /* start: 1-based */
static SegPos g_pos[MAX_SEG];

static int g_pe = 0;
static int g_klen = 0;          /* 拼接 tag 长度 */

/* ---------------- barcode hash (open addressing) ---------------- */
typedef struct {
    char  *keys;    /* cap * klen */
    int   *vals;    /* sample idx, 或 -2 = ambiguous */
    size_t cap, mask;
    int    klen;
} BarMap;

static void barmap_init(BarMap *m, int klen, size_t expect)
{
    size_t cap = 1024;
    while (cap < expect * 4) cap <<= 1;
    m->keys = xcalloc(cap, (size_t)klen);
    m->vals = xcalloc(cap, sizeof(int));
    m->cap = cap; m->mask = cap - 1; m->klen = klen;
}
/* 返回: 0=新插入 1=已存在相同 -1=冲突(不同样本) */
static int barmap_put(BarMap *m, const char *k, int v)
{
    size_t i = (size_t)(fnv(k, (size_t)m->klen) & m->mask);
    for (;;) {
        char *slot = m->keys + (size_t)i * m->klen;
        if (slot[0] == 0) {
            memcpy(slot, k, (size_t)m->klen);
            m->vals[i] = v;
            return 0;
        }
        if (memcmp(slot, k, (size_t)m->klen) == 0) {
            if (m->vals[i] == v) return 1;
            m->vals[i] = -2;      /* ambiguous */
            return -1;
        }
        i = (i + 1) & m->mask;
    }
}
static int barmap_get(const BarMap *m, const char *k)
{
    size_t i = (size_t)(fnv(k, (size_t)m->klen) & m->mask);
    for (;;) {
        const char *slot = m->keys + (size_t)i * m->klen;
        if (slot[0] == 0) return -1;
        if (memcmp(slot, k, (size_t)m->klen) == 0) return m->vals[i];
        i = (i + 1) & m->mask;
    }
}

static BarMap g_nei, g_exact;
static long   g_ambig_warn = 0;

/* ---------------- 邻居序列展开 (≤e mismatches, 替换为 != 原碱基) ---------------- */
typedef struct { char **v; int n, cap; } StrVec;

static void strvec_push(StrVec *sv, const char *s, int len)
{
    if (sv->n == sv->cap) {
        sv->cap = sv->cap ? sv->cap * 2 : 64;
        sv->v = xrealloc(sv->v, (size_t)sv->cap * sizeof(char *));
    }
    char *d = xmalloc((size_t)len + 1);
    memcpy(d, s, (size_t)len); d[len] = 0;
    sv->v[sv->n++] = d;
}
static void gen_rec(const char *seg, int len, int pos, int e, char *tmp, StrVec *out)
{
    if (pos == len) { strvec_push(out, tmp, len); return; }
    if (e > 0) {
        char o = seg[pos];
        for (const char *b = "ACGTN"; *b; b++) {
            if (*b == o) continue;
            tmp[pos] = *b;
            gen_rec(seg, len, pos + 1, e - 1, tmp, out);
        }
        tmp[pos] = o;
    }
    gen_rec(seg, len, pos + 1, e, tmp, out);
}
static void gen_variants(const char *seg, int len, int e, StrVec *out)
{
    char *tmp = xmalloc((size_t)len + 1);
    memcpy(tmp, seg, (size_t)len); tmp[len] = 0;
    for (int d = 0; d <= e; d++)
        gen_rec(seg, len, 0, d, tmp, out);
    free(tmp);
}

/* ---------------- 读入 barcode 列表 ---------------- */
static void load_barcodes(void)
{
    FILE *f = fopen(opt_bl, "r");
    if (!f) die("cannot open barcode list %s: %s", opt_bl, strerror(errno));

    int cap = 256;
    g_bc = xmalloc((size_t)cap * sizeof(Barcode));
    char line[4096];
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == 0) continue;
        char *fld[1 + MAX_SEG + 1];
        int nf = 0;
        char *tok = strtok(p, " \t\r\n");
        while (tok && nf < MAX_SEG + 1) { fld[nf++] = tok; tok = strtok(NULL, " \t\r\n"); }
        if (nf < 2) continue;
        int ns = nf - 1;
        if (g_nseg == 0) g_nseg = ns;
        if (ns != g_nseg)
            die("barcode list 第 %d 行序列列数(%d)与首行(%d)不一致", g_nbc + 1, ns, g_nseg);
        if (ns > MAX_SEG) die("barcode 段数最多支持 %d", MAX_SEG);
        if (g_nbc == cap) { cap *= 2; g_bc = xrealloc(g_bc, (size_t)cap * sizeof(Barcode)); }
        Barcode *b = &g_bc[g_nbc];
        snprintf(b->name, MAX_NAME, "%s", fld[0]);
        for (int i = 0; i < ns; i++) {
            /* mgi-dual: 交换两列, 使 tag/key 顺序与 Perl dual(第2列在前)一致 */
            char *s = (opt_mgi && ns == 2) ? fld[2 - i] : fld[i + 1];
            int L = (int)strlen(s);
            if (L == 0 || L >= MAX_SEQLEN) die("barcode 序列长度非法: %s", s);
            for (int j = 0; j < L; j++) s[j] = (char)toupper((unsigned char)s[j]);
            memcpy(b->seg[i], s, (size_t)L);
            b->seg[i][L] = 0;
            b->slen[i] = L;
        }
        g_nbc++;
    }
    fclose(f);
    if (g_nbc == 0) die("barcode list 为空: %s", opt_bl);

    if (opt_rc) {
        for (int i = 0; i < g_nbc; i++)
            for (int j = 0; j < g_nseg; j++)
                revcomp(g_bc[i].seg[j], g_bc[i].slen[j]);
    }
    g_klen = 0;
    for (int j = 0; j < g_nseg; j++) g_klen += g_bc[0].slen[j];
}

/* ---------------- 解析 barcode 位置 ---------------- */
static void parse_positions(void)
{
    for (int j = 0; j < g_nseg; j++) g_pos[j].len = g_bc[0].slen[j];

    if (opt_pos[0]) {
        char *spec = xstrdup(opt_pos);
        int idx = 0;
        char *tok = strtok(spec, ",");
        while (tok) {
            if (idx >= g_nseg) die("-p 给出的位置段数多于 barcode 列数(%d)", g_nseg);
            char *colon = strchr(tok, ':');
            if (colon) {
                *colon = 0;
                g_pos[idx].len = atoi(colon + 1);
                if (g_pos[idx].len <= 0) die("-p 长度非法: %s", tok);
            }
            g_pos[idx].start = atoi(tok);
            if (g_pos[idx].start <= 0) die("-p 位置非法(须为 1-based 正整数): %s", tok);
            idx++;
            tok = strtok(NULL, ",");
        }
        free(spec);
        if (idx != g_nseg) die("-p 给出的位置段数(%d)与 barcode 列数(%d)不一致", idx, g_nseg);
    } else {
        if (opt_first <= 0) die("必须提供 -f/--firstCycle 或 -p/--positions");
        if (opt_mgi && g_nseg == 2) {
            /* MGI_SplitDualBarcodes.pl 兼容: 列已交换(seg0=列表第2列), seg0 在 f, seg1 在 f+interval+1 */
            if (opt_interval < 0) die("双 barcode 必须提供 -i/--interval");
            g_pos[0].start = opt_first;
            g_pos[1].start = opt_first + opt_interval + 1;
        } else {
            g_pos[0].start = opt_first;
            if (g_nseg >= 2) {
                if (opt_interval < 0) die("双 barcode 必须提供 -i/--interval 或使用 -p 指定位置");
                g_pos[1].start = opt_first + g_pos[0].len + opt_interval;
            }
        if (g_nseg > 2) die("多于 2 段 barcode 必须使用 -p 指定各段位置");
    }
    /* -l/--lastCycle: MGI_SplitBarcode.pl 风格, 单 barcode 用 f..l 指定区间 */
    if (opt_last > 0) {
        if (g_nseg != 1) die("-l/--lastCycle 仅用于单 barcode 模式");
        if (opt_pos[0]) die("-l 与 -p 不可同时使用");
        if (opt_first <= 0) die("-l 需同时提供 -f");
        if (opt_last < opt_first) die("-l 必须 >= -f");
        g_pos[0].start = opt_first;
        g_pos[0].len = opt_last - opt_first + 1;
        if (g_pos[0].len != g_bc[0].slen[0])
            fprintf(stderr, "[split_barcode] WARNING: -f..-l 区间长度(%d)与列表 barcode 长度(%d)不一致, 将全部无法匹配\n",
                    g_pos[0].len, g_bc[0].slen[0]);
    }
    }
    /* 位置按 start 排序无关紧要, 但拼接 tag 按列表列顺序; 记录跨度用于切除 */
}

/* ---------------- chunk / 队列 ---------------- */
typedef struct Chunk {
    char  *arena;        /* 行内容连续存储 */
    size_t acap, aused;
    int   *loff, *llen;  /* 每行偏移/长度 */
    int    lcap, nlines;
    struct Chunk *next;
} Chunk;

typedef struct {
    Chunk *head, *tail;
    pthread_mutex_t mu;
    pthread_cond_t  cv;
    int  closed;
    long count;
} CQueue;

static CQueue g_work, g_pool;
static long   g_created = 0, g_pooled = 0;

static void cq_init(CQueue *q) { q->head = q->tail = NULL; q->closed = 0; q->count = 0;
    pthread_mutex_init(&q->mu, NULL); pthread_cond_init(&q->cv, NULL); }
static void cq_push(CQueue *q, Chunk *c)
{
    c->next = NULL;
    pthread_mutex_lock(&q->mu);
    if (q->tail) q->tail->next = c; else q->head = c;
    q->tail = c; q->count++;
    pthread_cond_signal(&q->cv);
    pthread_mutex_unlock(&q->mu);
}
static Chunk *cq_pop(CQueue *q)
{
    pthread_mutex_lock(&q->mu);
    while (!q->head && !q->closed) pthread_cond_wait(&q->cv, &q->mu);
    Chunk *c = q->head;
    if (c) { q->head = c->next; if (!q->head) q->tail = NULL; q->count--; }
    pthread_mutex_unlock(&q->mu);
    return c;
}
static void cq_close(CQueue *q)
{
    pthread_mutex_lock(&q->mu);
    q->closed = 1;
    pthread_cond_broadcast(&q->cv);
    pthread_mutex_unlock(&q->mu);
}

static void chunk_reset(Chunk *c) { c->aused = 0; c->nlines = 0; }

static int chunk_add_line(Chunk *c, const char *s, int len)
{
    if (c->aused + (size_t)len + 1 > c->acap) {
        while (c->aused + (size_t)len + 1 > c->acap) c->acap = c->acap ? c->acap * 2 : (1u << 22);
        c->arena = xrealloc(c->arena, c->acap);
    }
    if (c->nlines == c->lcap) {
        c->lcap = c->lcap ? c->lcap * 2 : 256;
        c->loff = xrealloc(c->loff, (size_t)c->lcap * sizeof(int));
        c->llen = xrealloc(c->llen, (size_t)c->lcap * sizeof(int));
    }
    int off = (int)c->aused;
    memcpy(c->arena + c->aused, s, (size_t)len);
    c->aused += (size_t)len;
    c->arena[c->aused++] = 0;
    c->loff[c->nlines] = off;
    c->llen[c->nlines] = len;
    return c->nlines++;
}

static Chunk *get_chunk(void)
{
    pthread_mutex_lock(&g_pool.mu);
    if (!g_pool.head && (g_created - g_pooled) >= max_inflight)
        while (!g_pool.head && (g_created - g_pooled) >= max_inflight)
            pthread_cond_wait(&g_pool.cv, &g_pool.mu);
    Chunk *c = g_pool.head;
    if (c) { g_pool.head = c->next; if (!g_pool.head) g_pool.tail = NULL; g_pooled--; }
    else { g_created++; c = xcalloc(1, sizeof(Chunk)); }
    pthread_mutex_unlock(&g_pool.mu);
    chunk_reset(c);
    return c;
}
static void put_chunk(Chunk *c)
{
    pthread_mutex_lock(&g_pool.mu);
    c->next = g_pool.tail ? g_pool.tail->next : NULL;
    if (g_pool.tail) { g_pool.tail->next = c; g_pool.tail = c; }
    else { g_pool.head = g_pool.tail = c; }
    g_pooled++;
    pthread_cond_signal(&g_pool.cv);
    pthread_mutex_unlock(&g_pool.mu);
}

/* ---------------- gz 逐行读取 ---------------- */
static int read_line(gzFile f, char **buf, size_t *cap)
{
    if (!*buf) { *cap = 8192; *buf = xmalloc(*cap); }
    size_t len = 0;
    for (;;) {
        if (gzgets(f, *buf + len, (int)(*cap - len)) == NULL) {
            if (len == 0) return -1;
            break;
        }
        len += strlen(*buf + len);
        if (len > 0 && (*buf)[len - 1] == '\n') break;
        if (len + 1 >= *cap) { *cap *= 2; *buf = xrealloc(*buf, *cap); }
        else break;
    }
    while (len > 0 && ((*buf)[len - 1] == '\n' || (*buf)[len - 1] == '\r')) len--;
    (*buf)[len] = 0;
    return (int)len;
}

/* ---------------- 样本输出(每文件一个写线程) ---------------- */
#define BLOCK_CAP (1u << 20)

typedef struct SBlock {
    struct SBlock *next;
    size_t len;
    char  *data;
} SBlock;

typedef struct {
    gzFile  f;
    SBlock *cur;            /* 正在填充的块 */
    SBlock *qhead, *qtail;  /* 待压缩块队列 */
    SBlock *pool;           /* 空闲块池 */
    pthread_mutex_t mu;
    pthread_cond_t  cv;
    int     eof;
    pthread_t tid;
} Sink;

typedef struct {
    char   name[MAX_NAME];
    Sink   s1, s2;
    int    dual;
} Sample;

static Sample g_unbar;                /* 未匹配 reads */
static Sample *g_sample = NULL;       /* 每个 barcode 一个 */

static SBlock *block_get(Sink *s)
{
    if (s->pool) {
        SBlock *b = s->pool;
        s->pool = b->next;
        b->len = 0;
        return b;
    }
    SBlock *b = xmalloc(sizeof(SBlock));
    b->data = xmalloc(BLOCK_CAP);
    b->len = 0;
    b->next = NULL;
    return b;
}

static void *sink_thread(void *arg)
{
    Sink *s = (Sink *)arg;
    for (;;) {
        pthread_mutex_lock(&s->mu);
        while (!s->qhead && !s->eof) pthread_cond_wait(&s->cv, &s->mu);
        SBlock *b = s->qhead;
        if (b) { s->qhead = b->next; if (!s->qhead) s->qtail = NULL; }
        pthread_mutex_unlock(&s->mu);
        if (!b) break;                     /* eof 且队列已空 */
        if (gzwrite(s->f, b->data, (unsigned)b->len) <= 0)
            fprintf(stderr, "[split_barcode] WARNING: gzwrite 失败\n");
        pthread_mutex_lock(&s->mu);
        b->next = s->pool; s->pool = b;
        pthread_mutex_unlock(&s->mu);
    }
    return NULL;
}

static void sink_write(Sink *s, const char *d, size_t n)
{
    pthread_mutex_lock(&s->mu);
    if (!s->cur) s->cur = block_get(s);
    if (s->cur->len + n > BLOCK_CAP) {
        /* 满块入队, 取新块 */
        SBlock *b = s->cur;
        if (s->qtail) s->qtail->next = b; else s->qhead = b;
        s->qtail = b; b->next = NULL;
        pthread_cond_signal(&s->cv);
        s->cur = block_get(s);
    }
    if (n > BLOCK_CAP) {
        /* 超大记录: 独立块直接入队(罕见路径) */
        SBlock *b = xmalloc(sizeof(SBlock));
        b->data = xmalloc(n);
        b->len = n;
        memcpy(b->data, d, n);
        b->next = NULL;
        if (s->qtail) s->qtail->next = b; else s->qhead = b;
        s->qtail = b;
        pthread_cond_signal(&s->cv);
    } else {
        memcpy(s->cur->data + s->cur->len, d, n);
        s->cur->len += n;
    }
    pthread_mutex_unlock(&s->mu);
}

static void sink_open(Sink *s, const char *path)
{
    memset(s, 0, sizeof(Sink));
    char mode[8];
    /* opt_compress=N 时用 zlib 透明模式(wT), 写出无 gzip 头的纯文本 */
    snprintf(mode, sizeof(mode), "wb%d%s", opt_level, opt_compress ? "" : "T");
    s->f = gzopen(path, mode);
    if (!s->f) die("cannot open output %s", path);
    gzbuffer(s->f, 1 << 19);
    pthread_mutex_init(&s->mu, NULL);
    pthread_cond_init(&s->cv, NULL);
    if (pthread_create(&s->tid, NULL, sink_thread, s) != 0)
        die("cannot create writer thread");
}

static void sink_close(Sink *s)
{
    pthread_mutex_lock(&s->mu);
    if (s->cur && s->cur->len) {
        if (s->qtail) s->qtail->next = s->cur; else s->qhead = s->cur;
        s->qtail = s->cur; s->cur->next = NULL;
        s->cur = NULL;
    } else if (s->cur) {
        s->cur->next = s->pool; s->pool = s->cur; s->cur = NULL;
    }
    s->eof = 1;
    pthread_cond_broadcast(&s->cv);
    pthread_mutex_unlock(&s->mu);
    pthread_join(s->tid, NULL);
    gzclose(s->f);
    s->f = NULL;
}

static void sample_open(Sample *s, const char *name, const char *p1, const char *p2, int dual)
{
    snprintf(s->name, MAX_NAME, "%s", name);
    s->dual = dual;
    sink_open(&s->s1, p1);
    if (dual) sink_open(&s->s2, p2);
}
static void sample_close(Sample *s)
{
    sink_close(&s->s1);
    if (s->dual) sink_close(&s->s2);
}

/* ---------------- 工作线程上下文 ---------------- */
typedef struct {
    char **keys;
    uint64_t *cnt;
    size_t cap, mask, n;
} TagMap;

static void tagmap_init(TagMap *m, size_t cap0)
{
    size_t cap = 1024;
    while (cap < cap0) cap <<= 1;
    m->keys = xcalloc(cap, sizeof(char *));
    m->cnt  = xcalloc(cap, sizeof(uint64_t));
    m->cap = cap; m->mask = cap - 1; m->n = 0;
}
static void tagmap_grow(TagMap *m)
{
    size_t ncap = m->cap << 1;
    char **nk = xcalloc(ncap, sizeof(char *));
    uint64_t *nc = xcalloc(ncap, sizeof(uint64_t));
    size_t nmask = ncap - 1;
    for (size_t i = 0; i < m->cap; i++) {
        if (!m->keys[i]) continue;
        size_t j = (size_t)(fnv(m->keys[i], strlen(m->keys[i])) & nmask);
        while (nk[j]) j = (j + 1) & nmask;
        nk[j] = m->keys[i]; nc[j] = m->cnt[i];
    }
    free(m->keys); free(m->cnt);
    m->keys = nk; m->cnt = nc; m->cap = ncap; m->mask = nmask;
}
static void tagmap_add(TagMap *m, const char *k, uint64_t v)
{
    if ((m->n + 1) * 10 >= m->cap * 7) tagmap_grow(m);
    size_t i = (size_t)(fnv(k, strlen(k)) & m->mask);
    for (;;) {
        if (!m->keys[i]) { m->keys[i] = xstrdup(k); m->cnt[i] = v; m->n++; return; }
        if (strcmp(m->keys[i], k) == 0) { m->cnt[i] += v; return; }
        i = (i + 1) & m->mask;
    }
}
static void tagmap_merge(TagMap *dst, const TagMap *src)
{
    for (size_t i = 0; i < src->cap; i++)
        if (src->keys[i]) tagmap_add(dst, src->keys[i], src->cnt[i]);
}

typedef struct {
    int      wid;
    TagMap   tags;
    long    *correct, *corrected;
    long     total;
    char    *out1, *out2;          /* 输出记录组装缓冲 */
    size_t   ocap1, ocap2;
    char    *tag;                  /* tag 组装缓冲 */
    char    *head;                 /* head 重写缓冲 */
    size_t   hcap;
} WorkerCtx;

static WorkerCtx *g_ctx;
static pthread_t *g_tids;
static pthread_t  g_rtid;

/* ---------------- 读线程 ---------------- */
static void *reader_thread(void *arg)
{
    (void)arg;
    gzFile f1 = gzopen(opt_r1, "rb");
    if (!f1) die("cannot open %s", opt_r1);
    gzbuffer(f1, 1 << 20);
    gzFile f2 = NULL;
    if (g_pe) {
        f2 = gzopen(opt_r2, "rb");
        if (!f2) die("cannot open %s", opt_r2);
        gzbuffer(f2, 1 << 20);
    }
    int lpr = g_pe ? 8 : 4;        /* lines per record */
    char *b1 = NULL, *b2 = NULL; size_t c1 = 0, c2 = 0;

    Chunk *c = get_chunk();
    for (;;) {
        int n = read_line(f1, &b1, &c1);
        if (n < 0) break;
        chunk_add_line(c, b1, n);
        for (int i = 1; i < lpr; i++) {
            int m = (i < 4) ? read_line(f1, &b1, &c1) : read_line(f2, &b2, &c2);
            if (m < 0) {
                fprintf(stderr, "[split_barcode] WARNING: 文件尾部记录不完整, 已丢弃\n");
                break;
            }
            chunk_add_line(c, i >= 4 ? b2 : b1, m);
        }
        if (c->nlines % lpr) break;   /* 不完整记录 */
        if (c->nlines >= chunk_records * lpr) {
            cq_push(&g_work, c);
            c = get_chunk();
        }
    }
    cq_push(&g_work, c);
    cq_close(&g_work);
    gzclose(f1);
    if (f2) gzclose(f2);
    free(b1); free(b2);
    return NULL;
}

/* ---------------- 工作线程 ---------------- */
/* 无换行追加 */
static void apend(char **buf, size_t *cap, size_t *u, const char *s, int n)
{
    if (*u + (size_t)n + 1 > *cap) {
        while (*u + (size_t)n + 1 > *cap) *cap = *cap ? *cap * 2 : (1u << 16);
        *buf = xrealloc(*buf, *cap);
    }
    memcpy(*buf + *u, s, (size_t)n);
    *u += (size_t)n;
}
/* 追加一行(带换行) */
static void oappend(char **buf, size_t *cap, size_t *u, const char *s, int n)
{
    apend(buf, cap, u, s, n);
    (*buf)[(*u)++] = '\n';
}

/* head 重写: 首 '/' 前部分 + "#UMI" + '/' + 后缀; 无 '/' 则直接追加 */
static void build_head(const char *h, int hlen, const char *umi, int ulen,
                       char **buf, size_t *cap, size_t *u)
{
    size_t need = (size_t)hlen + (size_t)ulen + 4;
    if (need > *cap) { *cap = need * 2; *buf = xrealloc(*buf, *cap); }
    char *p = *buf;
    const char *slash = memchr(h, '/', (size_t)hlen);
    size_t w = 0;
    if (slash) {
        size_t pre = (size_t)(slash - h);
        memcpy(p, h, pre); w = pre;
        p[w++] = '#';
        memcpy(p + w, umi, (size_t)ulen); w += (size_t)ulen;
        p[w++] = '/';
        memcpy(p + w, slash + 1, (size_t)(hlen - pre - 1));
        w += (size_t)(hlen - pre - 1);
    } else {
        memcpy(p, h, (size_t)hlen); w = (size_t)hlen;
        p[w++] = '#';
        memcpy(p + w, umi, (size_t)ulen); w += (size_t)ulen;
    }
    *u = w;
}

static void *worker_thread(void *arg)
{
    WorkerCtx *w = (WorkerCtx *)arg;
    int lpr = g_pe ? 8 : 4;
    int bseq = (g_pe ? (opt_bcread - 1) * 4 : 0) + 1;   /* barcode read 的 seq 行号 */
    int minstart = g_pos[0].start - 1, maxend = 0;
    for (int j = 0; j < g_nseg; j++) {
        int s0 = g_pos[j].start - 1;
        if (s0 < minstart) minstart = s0;
        if (s0 + g_pos[j].len > maxend) maxend = s0 + g_pos[j].len;
    }
    int nmates = g_pe ? 2 : 1;

    Chunk *c;
    while ((c = cq_pop(&g_work)) != NULL) {
        const char *L[8];
        for (int base = 0; base + lpr <= c->nlines; base += lpr) {
            for (int i = 0; i < lpr; i++)
                L[i] = c->arena + c->loff[base + i];
            w->total++;
            const char *seq_bc = L[bseq];
            int slen_bc = c->llen[base + bseq];

            /* 提取 tag(按列表列顺序拼接各段) */
            int ok = 1, off = 0;
            for (int j = 0; j < g_nseg; j++) {
                int s0 = g_pos[j].start - 1;
                if (g_pos[j].len != g_bc[0].slen[j] ||   /* 位置长度与列表不一致→无法匹配 */
                    s0 + g_pos[j].len > slen_bc) { ok = 0; break; }
                memcpy(w->tag + off, seq_bc + s0, (size_t)g_pos[j].len);
                off += g_pos[j].len;
            }
            w->tag[g_klen] = 0;
            int sidx = ok ? barmap_get(&g_nei, w->tag) : -1;
            if (sidx < 0 || sidx == -2) {
                /* 未匹配(长度不足/未知/歧义): 原样写入 unbarcoded(-x 时只计数) */
                if (!opt_no_unbar) {
                    size_t u = 0;
                    for (int i = 0; i < 4; i++)
                        oappend(&w->out1, &w->ocap1, &u, L[i], c->llen[base + i]);
                    sink_write(&g_unbar.s1, w->out1, u);
                    if (g_pe) {
                        u = 0;
                        for (int i = 4; i < 8; i++)
                            oappend(&w->out2, &w->ocap2, &u, L[i], c->llen[base + i]);
                        sink_write(&g_unbar.s2, w->out2, u);
                    }
                }
                tagmap_add(&w->tags, w->tag, 1);
                continue;
            }
            int exact = (barmap_get(&g_exact, w->tag) == sidx);
            if (exact) w->correct[sidx]++; else w->corrected[sidx]++;
            tagmap_add(&w->tags, w->tag, 1);   /* TagStat 统计所有 reads(含已匹配), 对齐 Perl */

            /* UMI/头后缀: mgi-dual 模拟 Perl dual: '#' + read2 @ (fc+len2-1) 起 9 碱基;
             * 其余模式: barcode 区之后 maxend 起的 opt_umi 个碱基 */
            const char *umi;
            int ulen;
            if (opt_mgi && g_nseg == 2 && g_pe) {
                int off = opt_first + g_bc[0].slen[0] - 1;   /* Perl: fc+len2-1, len2=列表第2列 */
                if (off < 0) off = 0;
                if (off >= slen_bc) { umi = seq_bc; ulen = 0; }
                else {
                    umi = seq_bc + off;
                    ulen = (off + 9 <= slen_bc) ? 9 : slen_bc - off;
                }
            } else {
                umi = seq_bc + maxend;
                ulen = opt_umi;
                if (maxend >= slen_bc) ulen = 0;
                else if (maxend + ulen > slen_bc) ulen = slen_bc - maxend;
            }

            for (int m = 1; m <= nmates; m++) {
                int h = (m - 1) * 4;
                char **obuf  = (m == 1) ? &w->out1 : &w->out2;
                size_t *ocap = (m == 1) ? &w->ocap1 : &w->ocap2;
                size_t u = 0;

                if (ulen > 0) {
                    size_t hu = 0;
                    build_head(L[h], c->llen[base + h], umi, ulen,
                               &w->head, &w->hcap, &hu);
                    oappend(obuf, ocap, &u, w->head, (int)hu);
                } else {
                    oappend(obuf, ocap, &u, L[h], c->llen[base + h]);
                }
                if (m == opt_bcread) {
                    if (opt_mgi && g_nseg == 2) {
                        /* Perl dual: read2 截断到 barcode 之前 (substr(seq2,0,fc-1)) */
                        apend(obuf, ocap, &u, seq_bc, minstart);
                        (*obuf)[u++] = '\n';
                        oappend(obuf, ocap, &u, L[h + 2], c->llen[base + h + 2]);
                        const char *qual = L[h + 3];
                        apend(obuf, ocap, &u, qual, minstart);
                        (*obuf)[u++] = '\n';
                    } else if (opt_cut_before) {
                        /* 仅保留 barcode 之前部分 (MGI_SplitBarcode.pl 兼容) */
                        apend(obuf, ocap, &u, seq_bc, minstart);
                        (*obuf)[u++] = '\n';
                        oappend(obuf, ocap, &u, L[h + 2], c->llen[base + h + 2]);
                        const char *qual = L[h + 3];
                        apend(obuf, ocap, &u, qual, minstart);
                        (*obuf)[u++] = '\n';
                    } else {
                        /* 切除 barcode 区间 [minstart, maxend) */
                        apend(obuf, ocap, &u, seq_bc, minstart);
                        apend(obuf, ocap, &u, seq_bc + maxend, slen_bc - maxend);
                        (*obuf)[u++] = '\n';
                        oappend(obuf, ocap, &u, L[h + 2], c->llen[base + h + 2]);
                        const char *qual = L[h + 3];
                        int qlen = c->llen[base + h + 3];
                        apend(obuf, ocap, &u, qual, minstart);
                        apend(obuf, ocap, &u, qual + maxend, qlen - maxend);
                        (*obuf)[u++] = '\n';
                    }
                } else {
                    oappend(obuf, ocap, &u, L[h + 1], c->llen[base + h + 1]);
                    oappend(obuf, ocap, &u, L[h + 2], c->llen[base + h + 2]);
                    oappend(obuf, ocap, &u, L[h + 3], c->llen[base + h + 3]);
                }
                sink_write(m == 1 ? &g_sample[sidx].s1 : &g_sample[sidx].s2, *obuf, u);
            }
        }
        put_chunk(c);
    }
    return NULL;
}

/* ---------------- 统计输出 ---------------- */
typedef struct { const char *key; uint64_t cnt; } TagRow;
static int tagrow_cmp(const void *a, const void *b)
{
    const TagRow *x = a, *y = b;
    if (x->cnt != y->cnt) return (x->cnt < y->cnt) ? 1 : -1;
    return strcmp(x->key, y->key);
}

static void write_stats(long total_reads)
{
    char path[4096];

    /* BarcodeStat.txt */
    snprintf(path, sizeof(path), "%s/BarcodeStat.txt", opt_outdir);
    FILE *bs = fopen(path, "w");
    if (!bs) die("cannot write %s", path);
    fprintf(bs, "#SpeciesNO\tCorrect\tCorrected\tTotal\tPct\n");
    long tc = 0, tco = 0, tt = 0;
    double tp = 0;
    for (int i = 0; i < g_nbc; i++) {
        long c = 0, co = 0;
        for (int w = 0; w < opt_threads; w++) { c += g_ctx[w].correct[i]; co += g_ctx[w].corrected[i]; }
        long tot = c + co;
        double pct = total_reads ? (double)tot / (double)total_reads * 100.0 : 0.0;
        fprintf(bs, "%s\t%ld\t%ld\t%ld\t%.4f%%\n", g_bc[i].name, c, co, tot, pct);
        tc += c; tco += co; tt += tot; tp += pct;
    }
    if (opt_mgi)
        fprintf(bs, "Total\t%ld\t%ld\t%ld\t%.4f%%\t%ld\n", tc, tco, tt, tp, total_reads);
    else
        fprintf(bs, "Total\t%ld\t%ld\t%ld\t%.4f%%\n", tc, tco, tt, tp);
    fclose(bs);

    /* TagStat.txt */
    TagMap all;
    tagmap_init(&all, 1 << 16);
    for (int w = 0; w < opt_threads; w++) tagmap_merge(&all, &g_ctx[w].tags);

    TagRow *rows = xmalloc((all.n ? all.n : 1) * sizeof(TagRow));
    size_t nr = 0;
    for (size_t i = 0; i < all.cap; i++)
        if (all.keys[i]) { rows[nr].key = all.keys[i]; rows[nr].cnt = all.cnt[i]; nr++; }
    qsort(rows, nr, sizeof(TagRow), tagrow_cmp);

    snprintf(path, sizeof(path), "%s/TagStat.txt", opt_outdir);
    FILE *ss = fopen(path, "w");
    if (!ss) die("cannot write %s", path);
    fprintf(ss, "#Sequence\tSpeciesNO\treadCount\tPct\n");
    for (size_t i = 0; i < nr; i++) {
        int sidx = barmap_get(&g_nei, rows[i].key);
        const char *nm = (sidx >= 0 && sidx != -2) ? g_bc[sidx].name : "unknown";
        double pct = total_reads ? (double)rows[i].cnt / (double)total_reads * 100.0 : 0.0;
        fprintf(ss, "%s\t%s\t%llu\t%.2f%%\n", rows[i].key, nm,
                (unsigned long long)rows[i].cnt, pct);
    }
    fclose(ss);

    for (size_t i = 0; i < all.cap; i++) free(all.keys[i]);
    free(all.keys); free(all.cnt); free(rows);
}

/* ---------------- 输出文件名 ---------------- */
static void make_prefix(const char *path, char *out, size_t n)
{
    const char *b = strrchr(path, '/');
    b = b ? b + 1 : path;
    snprintf(out, n, "%s", b);
    static const char *sufs[] = { ".fq.gz", ".fastq.gz", ".fq", ".fastq", ".gz", NULL };
    for (int i = 0; sufs[i]; i++) {
        size_t L = strlen(sufs[i]);
        size_t len = strlen(out);
        if (len > L && strcasecmp(out + len - L, sufs[i]) == 0) { out[len - L] = 0; break; }
    }
    size_t len = strlen(out);
    if (len > 2 && (strcmp(out + len - 2, "_1") == 0 || strcmp(out + len - 2, "_2") == 0))
        out[len - 2] = 0;
}

static void mkdir_p(const char *path)
{
    char tmp[4096];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (len && tmp[len - 1] == '/') tmp[len - 1] = 0;
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') { *p = 0; mkdir(tmp, 0755); *p = '/'; }
    }
    mkdir(tmp, 0755);
}

int main(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if      (match_opt(a, "-r1", "--read1")       && i + 1 < argc) opt_r1 = argv[++i];
        else if (match_opt(a, "-r2", "--read2")       && i + 1 < argc) opt_r2 = argv[++i];
        else if (match_opt(a, "-b",  "--barcodeList") && i + 1 < argc) opt_bl = argv[++i];
        else if (match_opt(a, "-f",  "--firstCycle")  && i + 1 < argc) opt_first = atoi(argv[++i]);
        else if (match_opt(a, "-i",  "--interval")    && i + 1 < argc) opt_interval = atoi(argv[++i]);
        else if (match_opt(a, "-e",  "--errNum")      && i + 1 < argc) opt_errnum = atoi(argv[++i]);
        else if (match_opt(a, "-w",  "--bcRead")      && i + 1 < argc) opt_bcread = atoi(argv[++i]);
        else if (match_opt(a, "-u",  "--umiLen")      && i + 1 < argc) opt_umi = atoi(argv[++i]);
        else if (match_opt(a, "-p",  "--positions")   && i + 1 < argc) snprintf(opt_pos, sizeof(opt_pos), "%s", argv[++i]);
        else if (match_opt(a, "-rc", "--revcom")      && i + 1 < argc) { const char *v = argv[++i]; opt_rc = (v[0]=='Y'||v[0]=='y'); }
        else if (match_opt(a, "-c",  "--compress")    && i + 1 < argc) { const char *v = argv[++i]; opt_compress = (v[0]=='Y'||v[0]=='y'); }
        else if (match_opt(a, "-z",  "--level")       && i + 1 < argc) opt_level = atoi(argv[++i]);
        else if (match_opt(a, "-l",  "--lastCycle")   && i + 1 < argc) opt_last = atoi(argv[++i]);
        else if (match_opt(a, "-x",  "--no-unbarcoded")) opt_no_unbar = 1;
        else if (match_opt(a, NULL,  "--cut-before"))    opt_cut_before = 1;
        else if (match_opt(a, "-t",  "--threads")     && i + 1 < argc) opt_threads = atoi(argv[++i]);
        else if (match_opt(a, NULL,  "--mgi-dual")) opt_mgi = 1;
        else if (match_opt(a, "-o",  "--outdir")      && i + 1 < argc) opt_outdir = argv[++i];
        else if (match_opt(a, "-h",  "--help")) { fprintf(stderr, "%s", g_usage); return 0; }
        else die("未知参数: %s\n%s", a, g_usage);
    }
    if (!opt_r1 || !opt_bl) die("缺少 -r1 或 -b\n%s", g_usage);
    if (opt_errnum < 0) opt_errnum = opt_r2 ? 2 : 1;
    g_pe = opt_r2 ? 1 : 0;
    if (!opt_bcread) opt_bcread = g_pe ? 2 : 1;
    if (opt_threads < 1) opt_threads = 1;
    if (opt_threads > 64) opt_threads = 64;
    if (max_inflight < opt_threads) max_inflight = opt_threads + 2;

    load_barcodes();
    parse_positions();

    /* 输出目录 */
    mkdir_p(opt_outdir);

    /* 构建 匹配哈希 */
    barmap_init(&g_exact, g_klen, (size_t)g_nbc);
    barmap_init(&g_nei,   g_klen, (size_t)g_nbc * 64);
    char *tag = xmalloc((size_t)g_klen + 1);
    for (int i = 0; i < g_nbc; i++) {
        int off = 0;
        for (int j = 0; j < g_nseg; j++) {
            memcpy(tag + off, g_bc[i].seg[j], (size_t)g_bc[i].slen[j]);
            off += g_bc[i].slen[j];
        }
        tag[g_klen] = 0;
        if (barmap_put(&g_exact, tag, i) != 0)
            die("barcode 重复: %s 与已有 barcode 完全相同", g_bc[i].name);
        barmap_put(&g_nei, tag, i);   /* 精确匹配 */
        for (int j = 0; j < g_nseg; j++) {
            StrVec sv = {0};
            gen_variants(g_bc[i].seg[j], g_bc[i].slen[j], opt_errnum, &sv);
            for (int v = 0; v < sv.n; v++) {
                if (strcmp(sv.v[v], g_bc[i].seg[j]) == 0) { free(sv.v[v]); continue; }
                int o = 0;
                for (int k = 0; k < g_nseg; k++) {
                    const char *src = (k == j) ? sv.v[v] : g_bc[i].seg[k];
                    memcpy(tag + o, src, (size_t)((k == j) ? (int)strlen(src) : g_bc[i].slen[k]));
                    o += (k == j) ? (int)strlen(src) : g_bc[i].slen[k];
                }
                tag[g_klen] = 0;
                if (barmap_put(&g_nei, tag, i) != 0) g_ambig_warn++;
                free(sv.v[v]);
            }
            free(sv.v);
        }
    }
    free(tag);
    if (g_ambig_warn)
        fprintf(stderr, "[split_barcode] WARNING: %ld 个容错邻居序列映射到多个 barcode, 这些序列按未匹配处理\n", g_ambig_warn);

    /* 打开输出 */
    char prefix[1024], p1[4096], p2[4096];
    const char *ext = opt_compress ? ".gz" : "";
    make_prefix(opt_r1, prefix, sizeof(prefix));
    g_sample = xcalloc((size_t)g_nbc, sizeof(Sample));
    for (int i = 0; i < g_nbc; i++) {
        if (g_pe) {
            snprintf(p1, sizeof(p1), "%s/%s_%s_1.fq%s", opt_outdir, prefix, g_bc[i].name, ext);
            snprintf(p2, sizeof(p2), "%s/%s_%s_2.fq%s", opt_outdir, prefix, g_bc[i].name, ext);
        } else {
            snprintf(p1, sizeof(p1), "%s/%s_%s.fq%s",  opt_outdir, prefix, g_bc[i].name, ext);
            p2[0] = 0;
        }
        sample_open(&g_sample[i], g_bc[i].name, p1, p2, g_pe);
    }
    if (g_pe) {
        snprintf(p1, sizeof(p1), "%s/%s_unbarcoded_1.fq%s", opt_outdir, prefix, ext);
        snprintf(p2, sizeof(p2), "%s/%s_unbarcoded_2.fq%s", opt_outdir, prefix, ext);
    } else {
        snprintf(p1, sizeof(p1), "%s/%s_unbarcoded.fq%s", opt_outdir, prefix, ext);
        p2[0] = 0;
    }
    if (!opt_no_unbar)
        sample_open(&g_unbar, "_unbarcoded", p1, p2, g_pe);

    fprintf(stderr, "[split_barcode] reads1=%s reads2=%s\n", opt_r1, opt_r2 ? opt_r2 : "-");
    fprintf(stderr, "[split_barcode] mode=%s barcodes=%d 段数=%d errNum=%d rc=%d threads=%d\n",
            g_pe ? "PE" : "SE", g_nbc, g_nseg, opt_errnum, opt_rc, opt_threads);

    /* 启动线程 */
    cq_init(&g_work); cq_init(&g_pool);
    g_ctx = xcalloc((size_t)opt_threads, sizeof(WorkerCtx));
    g_tids = xcalloc((size_t)opt_threads, sizeof(pthread_t));
    for (int w = 0; w < opt_threads; w++) {
        WorkerCtx *c = &g_ctx[w];
        c->wid = w;
        tagmap_init(&c->tags, 1 << 14);
        c->correct   = xcalloc((size_t)g_nbc, sizeof(long));
        c->corrected = xcalloc((size_t)g_nbc, sizeof(long));
        c->tag = xmalloc((size_t)g_klen + 1);
        if (pthread_create(&g_tids[w], NULL, worker_thread, c) != 0)
            die("cannot create worker thread");
    }
    if (pthread_create(&g_rtid, NULL, reader_thread, NULL) != 0)
        die("cannot create reader thread");

    pthread_join(g_rtid, NULL);
    for (int w = 0; w < opt_threads; w++) pthread_join(g_tids[w], NULL);

    long total = 0;
    for (int w = 0; w < opt_threads; w++) total += g_ctx[w].total;
    write_stats(total);

    for (int i = 0; i < g_nbc; i++) sample_close(&g_sample[i]);
    if (!opt_no_unbar) sample_close(&g_unbar);

    long matched = 0;
    for (int i = 0; i < g_nbc; i++)
        for (int w = 0; w < opt_threads; w++)
            matched += g_ctx[w].correct[i] + g_ctx[w].corrected[i];
    fprintf(stderr, "[split_barcode] total=%ld matched=%ld (%.2f%%) 输出目录: %s\n",
            total, matched, total ? (double)matched / total * 100.0 : 0.0, opt_outdir);
    return 0;
}

