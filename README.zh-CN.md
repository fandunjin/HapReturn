# HapReturn

**将合并挂载并调整后的组装，按原始来源重新分组。**

HapReturn 将原始来源组装比对到最终组装，利用覆盖证据恢复染色体来源。
支持 hifiasm hap1/hap2 输入，也支持任意数量、自定义名称的来源组装。
二倍体、四倍体、八倍体及更高倍性使用同一套流程。

## 安装与运行

需要 Python 3.9+；提供 FASTA 输入时需安装 minimap2。已有 PAF 可直接分析。

```bash
python -m pip install .
hapreturn --version
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt --threads 24 \
  --out-dir source_split
```

`chromosomes.txt` 每行一个染色体 ID。也可用 `--chromosome-count N` 选取前 N
条 FASTA 记录；均不指定时分析全部记录。未选入的序列保存在 `unplaced.fa`。
已有比对结果可改用 `--hap1-paf`、`--hap2-paf`；PAF 的 query 为原始来源 contig，
target 为最终组装。

## 八倍体及更高倍性

八倍体按两个来源回分时，可以追加：

```text
--ploidy 8 --homolog-groups groups.tsv --source-copies hap1=4 --source-copies hap2=4
```

`groups.tsv` 含 `group_id` 和 `seq_id` 两列，每行一条染色体；八倍体每个同源组
有 8 条记录。12 倍体、16 倍体可分别使用 `--ploidy 12`、`--ploidy 16` 并配置
对应的来源拷贝数。也支持组内来源数量不相等的情况。

提供多个已知来源时，使用制表符分隔的来源清单，例如：

```text
source  path         format
copy01  copy01.paf   paf
copy02  copy02.paf   paf
...
copy08  copy08.paf   paf
```

`format` 可选 `paf` 或 `fasta`，路径相对于清单文件。每个来源输出一个同名 FASTA。
运行内置八来源、两个同源组的八倍体示例：

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --sources examples/octoploid/sources.tsv \
  --homolog-groups examples/octoploid/groups.tsv --ploidy 8 \
  --min-aln-len 1 --out-dir octoploid_output
```

预期输出 8 套来源，每套 2 条染色体。1 bp 阈值仅用于此小型合成示例。
也可重复使用 `--source-paf 名称=路径` 或 `--source-fasta 名称=路径`。
软件按输入来源恢复分组：两个来源输入得到两套来源，多个来源输入得到对应多套。
倍性和来源数在分配逻辑中均无固定上限，名称由用户指定。

## 主要参数与输出

默认保留 mapQ ≥10、比对长度 ≥10 kb 的记录，合并重叠区间后计算来源覆盖量。
优势来源覆盖染色体至少 25%，且占来源覆盖量总和至少 60% 时接受该标签。
`--source-copies 名称=N` 检查每个同源组中的来源拷贝数，可重复指定。

输出包括各来源 FASTA、`unresolved.fa`、`unplaced.fa`、逐染色体证据表
`assignment.tsv`、同源组检查表 `group_qc.tsv` 和汇总表 `summary.tsv`。
保留染色体 ID，逐行写出序列以支持大型组装。

退出码 0 表示成功，1 表示输入或比对错误，2 表示有待审查染色体或组内数量不符；
退出码 2 仍保留输出。`--allow-unresolved` 接受部分结果，`--force` 允许复用非空
输出目录。软件按整条染色体分组，局部相位切换可在后续分析中检查。

完整参数见 [英文说明](README.md)，计算方法见 [方法说明](docs/method.md)。
采用 [MIT 许可证](LICENSE)。
