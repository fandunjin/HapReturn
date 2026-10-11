# HapReturn

**将合并调整后的组装，按原始 hifiasm hap1/hap2 支持与同源组约束回分为两套。**

二倍体、四倍体、八倍体及更高偶数倍性使用相同逻辑：每个同源组中两套各占一半，选择总来源支持最高的组合。

```text
original hap1 + original hap2 → merge → HapHiC / Juicebox → curated 2n
      |                |                                     |
      +----------------+------------ HapReturn --------------+
                                        |
                         alignments + source provenance
                                        |
                         homolog groups of P chromosomes
                                        |
                     maximize support with P/2 per source
                               /                \
                          hap1 (n)           hap2 (n)
                               \                /
                        PDF figures + HTML report
```

流程为：原始两套 contig 与调整后的序列比对 → 统计来源支持 → 每组 P 条染色体按总支持最高原则各分 P/2 条 → 输出两套序列和可视化报告，即 **2n → n + n**。

[English](README.md) · [方法说明](docs/method.md) · [图片示例](#图片示例)

## 主要功能

- **恢复两套序列：**利用原始 hifiasm hap1/hap2 的来源证据分配染色体。
- **组内等数分配：**满足每组各半的条件，选择总支持最高的组合。
- **追踪局部来源：**在连续灰色区段内扩大证据窗口，保留判定记录。
- **可视化与统计：**输出染色体数量、长度、分配表、三个独立 PDF 文件及一个完整 HTML 报告。

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

`window_qc.tsv` 列出扩增前后数量及剩余灰色的原因。扩大统计范围不会增加比对证据；整段缺少覆盖或两个来源竞争时，扩大到全长后仍可能未确定。更新这些图形可直接复用已有 PAF。

## 输出

| 文件 | 内容 |
|---|---|
| `hap1.fa`、`hap2.fa` | 两套分组序列，保留原始染色体 ID |
| `summary.tsv`、`assignment.tsv` | 数量长度统计、逐染色体证据与最终分配 |
| `homolog_groups.tsv`、`group_qc.tsv` | 实际同源组及约束分配检查 |
| `partition.pdf` | 来源支持与输出统计；染色体来源图 |
| `alignment_matrix.pdf` | 最终 hap1、hap2 各一张连续矩阵，随记录数调整尺寸 |
| `window_qc.tsv` | 局部窗口扩增前后统计及剩余灰色的原因 |
| `hap_comparison.pdf` | 每套 hap 一张连续对应关系图，随记录数增加图高 |
| `report.html` | 全部图形、来源统计、窗口悬停与编号筛选表 |

PDF 使用可编辑字体与独立文件；蓝、橙表示原始来源，灰色表示局部未确定证据。对应图保留两个原始来源的有效连接，每套最终 hap 连续展示，不再按记录数切块。

## 图片示例

以下图形全部使用**合成序列和模拟比对坐标**，展示三个四倍体同源组，最终两套各 6 条染色体。蓝色和橙色始终表示原始 hifiasm 来源，灰色表示局部未确定证据；整条染色体的最终分组与局部来源颜色分别展示。

### 来源支持与两套数量统计

逐染色体展示两个原始来源的独立覆盖比例，并统计输出序列的数量和长度；淡横线辅助阅读。这是 `partition.pdf` 的第一张图。

![合成数据：来源支持及恢复染色体数量](docs/images/source_support.png)

### 染色体来源图

按真实序列长度比例绘制染色体轮廓，hap1 和 hap2 分行展示，内部色带标记局部来源。示例中的 chr10 独立证据更支持 hap2，但在组内等数约束下分入 hap1；局部颜色仍保留原始支持信息。这是 `partition.pdf` 的第二张图。

![合成数据：染色体尺度的来源分布](docs/images/chromosome_provenance.png)

### 比对矩阵热图

行为原始 contig，列为调整后的染色体；单元格表示该 contig 对目标染色体的覆盖比例。最终 hap1、hap2 各自一张连续矩阵，保留两个原始来源的全部 contig 行，不按行列切页。侧边和顶部色标分别标明原始来源与最终归属，输出为 `alignment_matrix.pdf`；下图为 hap1 示例。

![合成数据：原始 contig 与调整后染色体的比对矩阵](docs/images/alignment_matrix.png)

### 原始与恢复后的对应关系

带状连线保留原始 contig 与调整后染色体的比对坐标，每套最终 hap 连续展示在一张图中，输出为 `hap_comparison.pdf`。下图为恢复后的 hap1；连线颜色表示原始来源，包含原始 hap2 对最终 hap1 染色体的有效比对。

![合成数据：原始 hifiasm 与恢复后 hap1 的对应关系](docs/images/hap_comparison.png)

### 复现图片示例

安装后，在软件仓库目录运行以下命令。脚本直接生成模拟 PAF，此示例无需 minimap2。

```bash
python examples/visualization/generate.py demo_input
hapreturn \
  --final-fasta demo_input/final.fa \
  --hap1-paf demo_input/hap1.paf --hap2-paf demo_input/hap2.paf \
  --chromosome-count 12 --ploidy 4 --plot-window-size 12000 \
  --allow-unresolved --out-dir demo_report
```

打开 `demo_report/report.html` 可查看全部图形和带编号的分配表。示例特意包含一个需复核的同源组，详见 `group_qc.tsv`；最终染色体数量为 6+6。

## 只更新报告

```bash
hapreturn-report --result-dir result --out-dir refreshed_report
```

此命令保留既有染色体分配，仅更新图形和局部窗口判定。要启用新的组内约束分配，请用原始最终 FASTA 和已有 PAF 重新运行主命令，输出到新目录。自定义过滤与窗口参数应与原分析一致。

## 独立绘制对应关系

```bash
hapreturn-chromosomes \
  --final-fasta curated.fa --hap1-paf hap1.paf --hap2-paf hap2.paf \
  --assignment result/assignment.tsv --out-pdf original_to_curated.pdf
```

退出码 0 表示通过，1 表示输入或运行错误，2 表示已有输出但存在需复核的证据或同分解。`--allow-unresolved` 可接受需复核输出。方法见 [Method](docs/method.md)，采用 MIT 许可证。
