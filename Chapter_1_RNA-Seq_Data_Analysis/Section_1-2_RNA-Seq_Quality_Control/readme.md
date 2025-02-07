**Chapter 1** RNA-Seq 数据处理与分析：比对、质量控制与定量（RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification）

**Section 1-2** RNA-Seq 质量控制（RNA-Seq quality control）

本节参考 *STAT 115 2021 Homework 1 Problem 2,3*

## 目录
  - [1 使用 RSeQC  进行质控](#1-使用-rseqc--进行质控)
    - [（1）查看 BAM 文件的基本统计信息](#1查看-bam-文件的基本统计信息)
    - [（2）RNA-seq 质控 (RNA-seq quality control)](#2rna-seq-质控-rna-seq-quality-control)
    - [（3）Gene Body Coverage 图](#3gene-body-coverage-图)
  - [2 用Python或R绘制 Gene Body Coverage 图](#2-用python或r绘制-gene-body-coverage-图)



我们在完成**比对**（alignment）后，接下来需要对比对结果进行**质量控制**（Quality Control）。

本节我们使用`RSeQC`工具进行质控，主要评估转录本（Transcript）的完整性和基因组覆盖情况。

## 1 使用 RSeQC  进行质控

需要用到的工具为`RSeQC`，RSeQC官网：[https://rseqc.sourceforge.net/](https://rseqc.sourceforge.net/)。
注意RSeQC不支持Windows，需要安装在Lunix系统上。

安装命令：

```bash
pip install RSeQC
```

或在conda环境下：

```bash
conda install RSeQC
```

![Image](https://github.com/user-attachments/assets/8b96bafd-fe3c-4c7b-b2d7-6686b95453d4)

安装好了之后，我们可以使用`RSeQC`下各个命令行工具了。RSeQC的各个命令行工具都是由Python的`py`文件构成。

### （1）查看 BAM 文件的基本统计信息

上一节中我们讲到了使用 `samtools flagstat` 和 `pysam`库 的方法查看比对结果`BAM`文件，现在我们也可以用RSeQC进行BAM统计信息的查看：

```bash
bam_stat.py -i 1M_SRR9336468_Aligned.sortedByCoord.out.bam
```

![Image](https://github.com/user-attachments/assets/fcf554ec-e8d9-42e0-b4d4-607561ae78c6)

可以看到`bam_stat.py`提供了更细颗粒度的统计信息。包括：**mapq < mapq_cut**：比对质量分数（MAPQ）低于阈值的非唯一比对 reads 数量；**mapq >= mapq_cut**：比对质量分数高于阈值的唯一比对 reads 数量。

### （2）RNA-seq 质控 (RNA-seq quality control)

**1 首先获取BED格式的gene model文件**

对于转录本（Transcript）的完整性评估，我们可用RSeQC进行**TIN** (transcript integrity number)。

这里需要用到基因模型（gene model）的BED文件。获取BED文件有很多种方法。我们可以使用Galaxy（[Galaxy网站](https://usegalaxy.org/)），将基因组的GTF文件转换为BED文件。

下图展示了将前面提供的 *Saccharomyces_cerevisiae.R64-1-1.107.gtf* 转换为BED12格式的文件。

![Image](https://github.com/user-attachments/assets/64494971-a95b-4b92-93ab-0cfec86a6856)

此外，也可以自己将 GTF文件 转化为 BED文件 ，BED文件本质上是简化版的基因组注释信息，可在网上自行搜索相关方法。

**2 生成BAM文件对应的`bai`文件**

使用`samtools`生成`bai`索引文件。

```bash
cd ~/STAR_results/
samtools index 1M_SRR9336468_Aligned.sortedByCoord.out.bam
samtools index 1M_SRR9336471_Aligned.sortedByCoord.out.bam
samtools index 1M_SRR9336474_Aligned.sortedByCoord.out.bam
```

**3 将上面转换好的`.bed12`的文件后缀改成`.bed`，使用`tin.py`跑程序。**

代码如下。`-i`参数下可以是一个目录，也可以是单个bam文件。具体参考文档：[tin.py](https://rseqc.sourceforge.net/#tin-py)

```bash
tin.py -i ~/STAR_results -r ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.bed
```

> **TIN（Transcript Integrity Number，转录本完整性数值）** 和 **MedTIN**（Median Transcript Integrity Number，转录本完整性数值的中位数） 是评估 RNA-seq 数据中转录本完整性的重要指标。它们用于衡量 RNA 样品在 RNA-seq 实验中的质量和降解程度，特别是在分析 RNA 的降解水平或样品完整性时非常有用。
>
> 根据RSeQC[原文档](https://rseqc.sourceforge.net/#tin-py)解释：**TIN** (transcript integrity number，转录本完整性指数) 与传统RIN（RNA Integrity Number，RNA完整性指数）相对应。但是RIN虽然是Sample（或**转录组**transcriptome） 层面最常用的评估RNA质量的指标，但是存在一定的缺陷，TIN则克服了这些缺陷。TIN为每个转录本计算一个分数（范围：0 ≤ TIN ≤ 100）。**MedTIN**（所有转录本 TIN 的中位数）则用于衡量 样本层面 的RNA完整性。

**4 结束后，每一个bam文件会生成两个结果文件。**

一个结尾是summary.txt，一个结尾是tin.xls，如下图。

![Image](https://github.com/user-attachments/assets/faff0914-cd43-4982-9a7f-49dbb4a43ecf)

---

我汇总了 *summary.txt* 中3个 BAM比对结果 的TIN评分：

| Bam_file                                    | TIN(mean) | TIN(median) | TIN(stdev) |
| ------------------------------------------- | --------- | ----------- | ---------- |
| 1M_SRR9336468_Aligned.sortedByCoord.out.bam | 85.85     | 91.49       | 15.18      |
| 1M_SRR9336471_Aligned.sortedByCoord.out.bam | 88.09     | 92.94       | 13.58      |
| 1M_SRR9336474_Aligned.sortedByCoord.out.bam | 84.50     | 90.59       | 16.05      |



### （3）Gene Body Coverage 图

**Gene Body Coverage 图** 是 RNA-seq 数据分析中用于评估测序数据均匀性的一个常用可视化工具。它展示了在基因体（gene body）的不同位置上测序覆盖度的分布情况，可用于评估测序数据质量。

使用RSeQC的 `geneBody_coverage.py` 方法可生成图。输入3个及以上的BAM文件，会生成一个折线图（lineGraph）和一个热图（heatmap）。

具体详见RSeQC文档：[genebody-coverage.py](https://rseqc.sourceforge.net/#genebody-coverage-py)

```bash
cd ~/STAR_results/
geneBody_coverage.py -r ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.bed -i 1M_SRR9336468_Aligned.sortedByCoord.out.bam,1M_SRR9336471_Aligned.sortedByCoord.out.bam,1M_SRR9336474_Aligned.sortedByCoord.out.bam  -o output
```

![Image](https://github.com/user-attachments/assets/3ae850c4-254f-447c-9fad-a8d8a3ea8a06)

结束后，会生成以下文件：

![Image](https://github.com/user-attachments/assets/b2fbc6ea-d885-4a13-8c55-6c3ef78816b6)

我们打开 *curves* 和 *heatmap* 的PDF文件查看：

![Image](https://github.com/user-attachments/assets/230798fb-226e-45a4-86c8-41f6c509bd26)

可以看到3个样本的比对质量总体都很不错，*SRR9336474* 和 *SRR9336471* 比 *SRR9336468* 更好一些。

## 2 用Python或R绘制 Gene Body Coverage 图

`geneBody_coverage.py` 生成的 `output.geneBodyCoverage.txt` 中包含了各 *gene body* 百分位上标准化后的reads数值。如下图。

![Image](https://github.com/user-attachments/assets/0153910e-0ed7-4cc0-947e-9958e41e4931)

我们可以使用 *Python* 及 *R* 绘制各样本的geneBodyCoverage图。代码如下：

**Python**

```python
import pandas as pd
import matplotlib.pyplot as plt

with open('output.geneBodyCoverage.txt') as f:
    lines = f.readlines()

percentiles = {}
for i in range(1,len(lines)):
    # skip i = 0, this is the header
    line_parts = lines[i].replace('\n','').split('\t')
    file = line_parts[0]
    percs = [float(x) for x in line_parts[1:]]
    percentiles[file] = percs

df = pd.DataFrame(percentiles)
x = range(1, 101)

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
plt.suptitle('Gene Body Coverage', size = 15, weight = 'bold')
for i in range(len(df.columns)):
    col = df.columns[i]
    print(col)
    axs[i].plot(x, df[col]/1000)
    axs[i].set_xlabel("Percentile of gene body (5' -> 3')")
    axs[i].set_ylabel("Read # (in thousands)")
    axs[i].set_title(
        'Library ' + col.split('_')[1]
    )

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
```

![Image](https://github.com/user-attachments/assets/af916265-3534-4811-bf9d-867d15f18534)

---

**R**

```R
# 在 Jupyter Notebook 的 R 环境中设置图形输出尺寸
options(repr.plot.width = 15, repr.plot.height = 5)

# 读取文件
lines <- readLines("output.geneBodyCoverage.txt")

# 提取数据
percentiles <- list()
for (i in 2:length(lines)) { # 从第2行开始（跳过header）
  line_parts <- unlist(strsplit(lines[i], "\t"))
  file <- line_parts[1]
  percs <- as.numeric(line_parts[-1])
  percentiles[[file]] <- percs
}

# 转换为数据框
df <- as.data.frame(percentiles)

# 定义 x 轴（1到100）
x <- 1:100

# 加载绘图包
library(ggplot2)

# 设置绘图窗口为多图排列（1行3列）
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

# 绘制图表
for (i in names(df)) {
  col <- i
  plot(
    x, df[[col]] / 1000,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "Percentile of gene body (5' -> 3')",
    ylab = "Read # (in thousands)",
    main = paste("Library", strsplit(col, "_")[[1]][2])
  )
}
```

![Image](https://github.com/user-attachments/assets/8ef6ae0f-bc35-4d24-ac93-9a55bee96548)
