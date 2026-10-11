# HapReturn

**将合并调整后的 2n 组装，按原始 hifiasm hap1/hap2 来源回分为两套 n。**

软件将原始 hap1、hap2 contig 分别比对到最终组装，根据覆盖证据恢复两套来源。
无论二倍体、四倍体、八倍体或更高倍性，逻辑始终是 **2n → n + n**。
这里 n 是总染色体数的一半，例如 2n=64，目标就是两套各 32 条染色体。

## 安装与运行

需要 Python 3.9+；提供 FASTA 时需安装 minimap2，已有 PAF 可直接分析。

```bash
python -m pip install .
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt --threads 24 \
  --out-dir result
```

不同倍性使用同一条命令，无需指定倍性。`chromosomes.txt` 每行一个最终组装的
染色体 ID；也可用 `--chromosome-count N` 选择前 N 条记录，均不指定时分析全部。
未选入的记录保存在 `unplaced.fa`。

已有比对结果可改用 `--hap1-paf`、`--hap2-paf`；query 为原始来源 contig，
target 为最终组装。默认输出 `hap1.fa`、`hap2.fa` 两套来源。

## 示例

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --hap1-paf examples/octoploid/hap1.paf \
  --hap2-paf examples/octoploid/hap2.paf \
  --min-aln-len 1 --out-dir example_result
```

预期：16 条输入染色体回分为 hap1、hap2 各 8 条。1 bp 阈值仅用于此合成示例。

## 可视化

运行后自动生成四类图及四面板总览，均提供 PNG、SVG、PDF：

- 染色体来源覆盖度：并列显示两个来源的独立覆盖比例。
- 分配散点图：颜色表示分配结果，点大小表示染色体长度。
- 染色体窗口来源图：展示沿染色体的局部来源证据。
- 输出统计图：分别展示两套来源、待确认及未选入记录的数量与总长度。

![分型结果示例](docs/images/overview.png)

打开 `report.html` 可离线查看图表，悬停查看窗口坐标及两来源覆盖量，
并按名称或来源筛选染色体表。灰色窗口表示证据未能确定来源。
`window_provenance.tsv` 保存逐窗口去重覆盖量，坐标为从 0 开始的左闭右开区间。

`--plot-window-size N` 设置窗口宽度（bp），默认按最长染色体约 100 个窗口
自动选择；`--no-plots` 可关闭可视化。窗口使用与整条染色体相同的覆盖度和
优势比例阈值进行局部计算。所有倍性均使用同一套两来源图表。

## 结果与参数

默认保留 mapQ ≥10、比对长度 ≥10 kb 的记录，合并重叠区间后计算两个来源的
覆盖量。优势来源覆盖染色体至少 25%，且占两来源覆盖量之和至少 60% 时接受标签。
阈值可通过参数调整。

输出包括 `hap1.fa`、`hap2.fa`、待确认的 `unresolved.fa`、未选入的 `unplaced.fa`、
逐染色体证据表 `assignment.tsv` 和汇总表 `summary.tsv`。保留染色体 ID，流式
写出序列。来源按比对证据判定，数量不通过强行分配来凑齐。

`--ploidy`、`--homolog-groups` 和组内来源数量参数仅用于可选检查，不改变拆分逻辑。
也可为两个来源自定义名称。退出码 0 表示成功，1 表示输入或比对错误，2 表示
有待确认结果；退出码 2 时仍保留输出。

完整参数见 [英文说明](README.md)，计算方法见 [方法说明](docs/method.md)。
采用 [MIT 许可证](LICENSE)。
