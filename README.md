# Sylph BIL4 四类双歧杆菌定量方法

这是一个可直接交给合作者使用的、非 Docker 的 Sylph 0.5.1 分析包。它针对宏基因组双端 clean FASTQ，定量四个预先确定的双歧杆菌参考分支：

- `BI`: *Bifidobacterium infantis*
- `BL`: *Bifidobacterium longum* subsp. *longum*
- `BS`: *Bifidobacterium longum* subsp. *suis/suillum* clade
- `BX`: *Bifidobacterium longum* subsp. X

包内已经包含四个参考 FASTA 和预构建的 `BIL4.syldb`，因此使用者不需要自行下载参考基因组或重新建库。只需要准备 clean paired-end FASTQ，并填写一个输入清单。

> 本包按当前目录中实际运行的工具整理，软件名称是 `Sylph`（不是 `spyla`），版本固定为 0.5.1。当前包仅测试过 Linux x86_64。

## 1. 安装环境

推荐使用 Miniforge/Mambaforge 的 `mamba`；没有 `mamba` 时可以把下面命令中的 `mamba` 换成 `conda`。

```bash
mamba env create -f environment.yml
mamba activate sylph-bil4
bash scripts/check_installation.sh
```

安装依赖只有：

- `sylph=0.5.1`
- `python=3.11`

Python 脚本只使用标准库，不需要 `pip install`、R、HUMAnN、MetaPhlAn 或 GNU Parallel。

## 2. 准备输入

输入必须是已经完成质控的双端 FASTQ；支持未压缩的 `.fastq/.fq` 和 gzip 压缩的 `.fastq.gz/.fq.gz`。每个样本一行，清单必须是制表符分隔的 TSV：

```text
sample_id<TAB>r1<TAB>r2
sample_001<TAB>/absolute/path/sample_001.clean_1.fastq.gz<TAB>/absolute/path/sample_001.clean_2.fastq.gz
sample_002<TAB>/absolute/path/sample_002.clean_1.fastq.gz<TAB>/absolute/path/sample_002.clean_2.fastq.gz
```

上面的 `<TAB>` 表示一个真实的制表符，不是四个普通字符。

可复制模板：

```bash
cp examples/input_manifest.example.tsv my_samples.tsv
# 用文本编辑器把 my_samples.tsv 中的示例路径替换为真实路径
```

说明：

- `sample_id` 必须唯一，建议只使用字母、数字、点、下划线和短横线。
- `r1` 和 `r2` 必须存在且大小大于 0；相对路径相对于清单所在目录解释。
- R1/R2 必须来自同一个样本，不能把不同样本或不同 lane 的文件误配成一对。
- 不要把原始 FASTQ、个人路径、运行日志或中间文件提交到 GitHub。

## 3. 运行

在本方法包根目录运行：

```bash
bash scripts/run_sylph_bil4.sh \
  --manifest my_samples.tsv \
  --output results/my_run \
  --threads 8 \
  --jobs 4
```

参数含义：

- `--threads`: Sylph profile 阶段使用的线程数。
- `--jobs`: 同时进行多少个样本的 sketch；每个 sketch 固定使用 1 个 Sylph 线程。总资源大致由 `jobs + threads` 决定，请根据服务器资源调整。
- `--output`: 结果目录；如果任务中断，可以用同一命令重新运行，已完成的 sketch 会复用。

如果 `sylph` 不在当前环境的 `PATH` 中，可以显式指定：

```bash
SYLPH_BIN=/path/to/sylph bash scripts/run_sylph_bil4.sh \
  --manifest my_samples.tsv \
  --output results/my_run
```

## 4. 输出文件

最重要的交付文件是：

```text
results/my_run/processed/BIL4_abundance_wide.csv
```

建议一并返回：

```text
results/my_run/qc/BIL4_summary.csv
results/my_run/processed/BIL4_long.csv
results/my_run/run_info.txt
```

文件用途：

- `BIL4_abundance_wide.csv`: 每行一个样本，四类双歧杆菌的主定量矩阵。
- `BIL4_long.csv`: 每个样本 × 每个分支一行，包含 ANI、coverage 和是否检出的审计表。
- `BIL4_summary.csv`: 样本数、缺失配对、检出数、阈值检查和数值范围等 QC。
- `raw/BIL4_sylph_raw.tsv`: Sylph 原始 profile 输出，通常不需要发送给下游分析人员。
- `run_info.txt`: 运行参数、软件版本和参考数据库信息。

主矩阵中的字段：

- `BI_abundance`, `BL_abundance`, `BS_abundance`, `BX_abundance`: Sylph `Taxonomic_abundance`；如果未达到该分支的 ANI 阈值，则记为 `0`。
- `BI_ANI`, `BL_ANI`, `BS_ANI`, `BX_ANI`: Sylph `Adjusted_ANI`；未得到有效 profile 时为 `NA`。
- `*_detected`: 是否达到该分支的阈值，`1` 为检出，`0` 为未检出。
- `BIL_total_abundance`: 四个分支 abundance 之和。
- `*_fraction_in_BIL`: 在四个目标分支总 abundance 中的比例；如果四个分支总和为 0，则为 `NA`。

## 5. 固定方法参数

当前版本固定使用：

| 项目 | 设置 |
|---|---|
| Sylph | 0.5.1 |
| genome sketch subsampling | `-c 200` |
| Sylph 默认 k-mer | `k=31` |
| profile 最低 ANI | `-m 97` |
| BI / BL / BS 后处理阈值 | Adjusted ANI ≥97 |
| BX 后处理阈值 | Adjusted ANI ≥99 |
| 输入 | paired-end clean FASTQ |

注意：profile 阶段统一使用 `-m 97`，BX 的更严格阈值是在结果整理时执行。因此不要直接把 raw profile 中的 BX 行当作最终 BX 检出结果，应该使用 `BIL4_abundance_wide.csv` 中的 `BX_detected` 和 `BX_abundance`。

## 6. 常见问题

### `sylph: command not found`

确认已经激活环境：

```bash
mamba activate sylph-bil4
sylph --version
```

应显示 `sylph 0.5.1`。

### `no clean paired FASTQ` 或 `missing FASTQ`

检查 TSV 是否为制表符分隔、路径是否写对、文件是否有读取权限。可以先运行：

```bash
python scripts/validate_inputs.py --manifest my_samples.tsv
```

### `duplicate sample_id`

一个清单中不能有重复的 `sample_id`。如果一个生物样本有多个测序 lane，应先在上游合并 clean reads，或给每个独立样本使用不同 ID。

### 服务器不是 Linux x86_64

本包当前只验证了 Linux x86_64。Windows 建议使用 WSL2；macOS 或 ARM 服务器需要另外验证 Sylph 的可用二进制，不要直接把本包结果与 Linux 结果混用而不做测试。

## 7. 参考文件与复现

`database/sylph_BIL4/reference_manifest.tsv` 记录了四个参考文件、reference ID、阈值和 MD5。包内 `BIL4.syldb` 是由这四个 FASTA 使用 Sylph 0.5.1、`-c 200` 构建的预计算数据库。

如需重新建库：

```bash
bash scripts/build_database.sh
```

重新建库前请确认四个 FASTA 的 MD5 与 manifest 一致。除非参考文件发生变化，不建议重新建库；直接使用包内数据库最容易保证不同使用者得到一致结果。

## 8. 给方法使用者的最短说明

可以直接把下面几项发给运行人员：

1. 安装 `environment.yml` 中的环境。
2. 将自己的 clean paired FASTQ 写入 `sample_id/r1/r2` 三列 TSV。
3. 运行 `bash scripts/run_sylph_bil4.sh --manifest ... --output ...`。
4. 把 `processed/BIL4_abundance_wide.csv`、`qc/BIL4_summary.csv` 和 `run_info.txt` 发回。

方法代码、参考 FASTA 和数据库应与本包一起作为一个版本保存；不要只单独传一个结果矩阵而不记录方法版本。
