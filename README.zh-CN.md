# HapReturn

**将合并并经 HapHiC/Juicebox 调整的组装，按原始 hifiasm hap1/hap2 来源重新分组。**

名称中的 Return 表示回溯到原始来源。软件输出两套来源 FASTA，以及可审查的
覆盖证据和歧义记录。染色体名称、数量和倍性均不写死。

## 安装与运行

需要 Python 3.9+。进入本软件目录安装：

```bash
python -m pip install .
hapreturn --version
```

提供原始 FASTA 时，需要额外安装 minimap2 并加入 PATH：

```bash
hapreturn \
  --final-fasta out_JBAT.FINAL.fa \
  --hap1-fasta original.hap1.fa --hap2-fasta original.hap2.fa \
  --chromosome-list chromosomes.txt --threads 24 \
  --out-dir source_split
```

`chromosomes.txt` 每行一个最终组装中的染色体 ID，其余序列输出到 `unplaced.fa`。
也可以用 `--chromosome-count N` 指定 FASTA 中前 N 条记录；不指定时分析全部记录。
每个来源也可以通过 `--hap1-paf` / `--hap2-paf` 提供已有 PAF，无需调用 minimap2。
PAF 的 target 必须是最终组装，query 必须是原始来源 contig。

## 判断与倍性

软件分别合并两个来源在染色体上的比对区间，避免重叠区间重复计数。默认保留
mapQ ≥10、比对长度 ≥10 kb 的记录，排除 secondary hits；优势来源覆盖染色体
至少 25%，且占两来源覆盖量之和至少 60% 时接受该来源标签。

`--ploidy P` 用于检查染色体数是否能被 P 整除；提供含 `group_id`、`seq_id`
两列的 `--homolog-groups groups.tsv` 后，还检查每个同源组是否有 P 条染色体。
`--hap1-per-group K` 可检查组内 hap1 数目，不会强制改写来源。

例如四倍体可检查 P=4、K=2；有依据的三倍体 1:2 来源构成可检查 P=3、K=1。
倍性参数不会自动推断染色体对应关系，也不会把两个来源分成更多生物学单倍型。

## 输出与状态

- `hap1.fa`、`hap2.fa`：满足阈值的两套来源染色体，保留原 ID。
- `unresolved.fa`：低覆盖或来源不明确的染色体。
- `unplaced.fa`：未选入染色体分析的记录。
- `assignment.tsv`：每条染色体的来源覆盖量、优势比例和判定。
- `group_qc.tsv`：可选的同源组及来源数量检查。

退出码 0 表示运行通过，1 表示输入或比对错误，2 表示需要审查歧义或同源组检查。
退出码 2 时仍保留可审查输出；`--allow-unresolved` 可以接受部分结果。
默认要求输出目录为空，重复使用非空目录需加 `--force`。

## 使用边界

该版本按整条染色体的优势来源分类。混合来源 scaffold 和局部 phase switch
可能被整体标签掩盖；A/B 亚基因组或 A1/A2/B1/B2 分配仍需共线性等进一步证据。
FASTA 需为未压缩文件；软件将组装读入内存，所需 RAM 随组装大小增加。

0.1.0 经过合成数据测试，尚未进行真实组装准确性评估和真实 minimap2 路径验证。
完整参数、示例和方法见 [英文首页](README.md) 和 [方法说明](docs/method.md)。
采用 MIT 许可证。
