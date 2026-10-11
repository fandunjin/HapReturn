# HapReturn

**将合并调整后的组装，按原始 hifiasm hap1/hap2 支持与同源组约束回分为两套。**

二倍体、四倍体、八倍体及更高偶数倍性使用相同逻辑：每个同源组中两套各占一半，选择总来源支持最高的组合。

## 安装

需要 Python 3.9+；FASTA 比对需 minimap2，已有 PAF 可直接复用。

```bash
python -m pip install .
hapreturn --version
```

## 输入顺序要求

在 Juicebox 中将主要染色体放在最前面，并按同源组连续排列。四倍体为第 1–4 条一组、第 5–8 条一组，以此类推；八倍体为每 8 条一组。组内两套各半的假设应符合输入设计。倍性 P 表示每组的同源染色体数，不固定亚基因组名称。

`--chromosome-count N` 选择 FASTA 的前 N 条记录，按照文件顺序分组，不按长度重排。其余记录输出到 `unplaced.fa`。不连续的同源组可用 `--homolog-groups groups.tsv` 显式指定，表头为 `group_id`、`seq_id`，须列出全部入选染色体。

## 运行：四倍体 56 条染色体

```bash
hapreturn \
  --final-fasta curated.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-count 56 --ploidy 4 --threads 24 \
  --out-dir result
```

默认采用组内约束分配：每 4 条选 2 条给 hap1，其余给 hap2，最终两套各 28 条。八倍体使用 `--ploidy 8`，每组为 4+4。

已有比对可用 `--hap1-paf` 和 `--hap2-paf` 替代 FASTA 来源输入，无需重新运行 minimap2。PAF query 为原始 contig，target 为调整后的组装。来源可用 `--source-paf LABEL=PATH` 自定义命名，始终提供两套。

## 分配方法与证据

默认保留 mapQ ≥10、比对长度 ≥10 kb 的记录。每条染色体分别计算两个来源比对区间的并集覆盖碱基数 C1、C2。每个 P 条染色体的组内，将 C1−C2 最大的 P/2 条分给 hap1，其余分给 hap2；这就是满足等数约束时总覆盖支持最大的解。

`assignment.tsv` 同时保留最终 `assignment`、独立证据 `evidence_assignment`、`group_id` 和最终来源覆盖比例。独立证据的默认阈值为覆盖至少 25%、来源占比至少 60%。`group_qc.tsv` 记录总支持、最优解与次优解差值，以及需要复核的证据变化或同分解。同分时按输入顺序稳定分配并标记复核；零证据不被描述为确定的来源证据。

如需沿用独立来源恢复，使用 `--assignment-mode provenance`；此模式允许不等数结果，倍性与组表仅用于可选检查。

## 局部灰色区段的窗口扩增

对连续 unresolved 灰色区段，以原窗口为基础逐级扩大支持统计范围，直到该连续灰色区段全长。使用已有 PAF 的覆盖证据，不重新运行比对，也不借用区段外已确定来源的窗口。达到覆盖和优势阈值后赋予局部来源，仍无充分证据则保留灰色。

`window_provenance.tsv` 保留原始窗口覆盖、原始归属、最终局部归属、判定方法及证据范围；`window_rescue.tsv` 记录每次扩增尝试。局部扩增独立于染色体的组内约束分配。

## 输出

| 文件 | 内容 |
|---|---|
| `hap1.fa`、`hap2.fa` | 两套分组序列，保留原始染色体 ID |
| `summary.tsv`、`assignment.tsv` | 数量长度统计、逐染色体证据与最终分配 |
| `homolog_groups.tsv`、`group_qc.tsv` | 实际同源组及约束分配检查 |
| `partition.pdf` | 来源支持与输出统计；染色体来源图 |
| `alignment_matrix.pdf` | 原始 contig × 最终染色体的比对矩阵 |
| `hap_comparison.pdf` | 每套 hap 一张连续对应关系图，随记录数增加图高 |
| `report.html` | 全部图形、来源统计、窗口悬停与编号筛选表 |

PDF 使用可编辑字体与独立文件；蓝、橙表示原始来源，灰色表示局部未确定证据。对应图保留两个原始来源的有效连接，每套最终 hap 连续展示，不再按记录数切块。

### 只更新报告

```bash
hapreturn-report --result-dir result --out-dir refreshed_report
```

此命令保留既有染色体分配，仅更新图形和局部窗口判定。要启用新的组内约束分配，请用原始最终 FASTA 和已有 PAF 重新运行主命令，输出到新目录。自定义过滤与窗口参数应与原分析一致。

### 独立绘制对应关系

```bash
hapreturn-chromosomes \
  --final-fasta curated.fa --hap1-paf hap1.paf --hap2-paf hap2.paf \
  --assignment result/assignment.tsv --out-pdf original_to_curated.pdf
```

退出码 0 表示通过，1 表示输入或运行错误，2 表示已有输出但存在需复核的证据或同分解。`--allow-unresolved` 可接受需复核输出。方法见 [Method](docs/method.md)，采用 MIT 许可证。
