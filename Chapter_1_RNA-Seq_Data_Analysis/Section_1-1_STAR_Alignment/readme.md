**Chapter 1** RNA-Seq 数据处理与分析：比对、质量控制与定量（RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification）

**Section 1-1** 用 STAR 进行 RNA-seq 测序比对（STAR alignment）

本节参考 *STAT 115 2021 Homework 1 Problem 1*

## 目录
  - [1 安装STAR](#1-安装star)
  - [2 数据准备](#2-数据准备)
    - [准备3类数据文件](#准备3类数据文件)
    - [预览基因组信息](#预览基因组信息)
  - [3 索引生成](#3-索引生成)
  - [4 RNA-Seq 数据比对](#4-rna-seq-数据比对)
    - [（1）获取酵母 RNA-Seq 原始数据](#1获取酵母-rna-seq-原始数据)
    - [（2）测试运行 STAR](#2测试运行-star)
  - [5 RNA-Seq 比对结果查看](#5-rna-seq-比对结果查看)
    - [（1）查看日志文件的统计数据](#1查看日志文件的统计数据)
    - [（2）查看BAM文件](#2查看bam文件)
      - [①使用Linux的Samtools工具](#①使用linux的samtools工具)
      - [②使用Python库 pysam](#②使用python库-pysam)


本章起我们对高通量测序（High-Throughput Sequencing, HTS）数据，即 **RNA-Seq** 数据进行处理与分析，包括比对、质量控制与定量。本节是测序比对（**Alignment**）的实战，我们使用常用的`STAR`比对方法，采用酵母基因组（Saccharomyces cerevisiae）作为研究对象。

## 1 安装STAR

我们通过 Conda 创建一个新环境 *star_env*，然后安装`STAR`。

```bash
conda create -n star_env -c bioconda star
conda activate star_env
```

验证安装：

```bash
STAR --version
```

## 2 数据准备

### 准备3类数据文件

我们需要准备3类数据文件：

基因组参考文件、注释文件、RNA-Seq 原始数据。

- **基因组参考文件**：

  基因组参考文件代表一个完整的 DNA 序列，RNA-Seq 比对时使用它作为参考。

  来源常为一些公共数据库，如**Ensembl**，可以提供例如人类基因组 `GRCh38`。文件格式为*FASTA*。

- **注释文件**：

  注释文件与基因组参考文件对应。例如**Ensembl**中人类基因组对应的 `GRCh38.gtf` 。**GTF 文件（Gene Transfer Format）** 是一种常用的基因注释文件格式，主要用于描述基因组上基因和转录本的具体坐标及其特征。

- **RNA-Seq 原始数据**：

  RNA-Seq 原始数据（常为*FASTQ*格式）一般是测序平台（如 Illumina）生成的读序文件（例如 `sample_1.fastq.gz`），用于和基因组参考文件比对，也是我们跑STAR的目的。

---

**下载参考基因组**：

**（1）**创建目录存放数据：

```bash
mkdir -p ~/STAR_data/genome
cd ~/STAR_data/genome
```

**（2）**下载*ensembl*上的酵母基因组和注释文件，并解压缩（该数据源介绍见下文）：

```
wget ftp://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-107/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz
gunzip *.gz
```

### 预览基因组信息

我们可以使用Python的`biopython`库对基因组的FASTA文件进行读取，查看基本信息。此外也可以使用其他软件，比如`IGV`可以可视化全基因组。

`biopython`库安装：

```bash
pip install biopython
```

预览基因组信息：

```python
from Bio import SeqIO

# 定义文件路径
file_path = "genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

# 读取并打印基本信息
with open(file_path, "r") as file:
    records = list(SeqIO.parse(file, "fasta"))
    
    # 打印序列的概览信息
    print(f"共有 {len(records)} 条序列")
    for i, record in enumerate(records[:5]):  # 查看前 5 条记录
        print(f"\n序列 {i + 1}:")
        print(f"ID: {record.id}")
        print(f"描述: {record.description}")
        print(f"序列长度: {len(record.seq)}")
        print(f"前 100 个碱基: {record.seq[:100]}")
```

![Image](https://github.com/user-attachments/assets/ca385422-587a-4b5e-b4bb-8169723d1222)

```python
# 提取每条序列的长度和 ID
seq_lengths = [len(record.seq) for record in records]
seq_ids = [record.id for record in records]

import matplotlib.pyplot as plt
# 绘制条形图
plt.figure(figsize=(10, 6))
plt.barh(seq_ids, seq_lengths, color='skyblue')
plt.xlabel("Sequence Length (bp)")
plt.ylabel("Sequence ID")
plt.title("Overview of Saccharomyces cerevisiae Genome Sequences")
plt.tight_layout()
plt.show()
```

![Image](https://github.com/user-attachments/assets/1f4276ef-78c1-4b3b-bc52-f4837d551790)



## 3 索引生成

**生成索引（indexing）** 这一步的主要目的是：

1. 加速 reads 映射到参考基因组的过程，优化比对性能
2. 支持跨剪切位点比对（spliced alignment）

我们使用参考基因组和注释文件，来生成 STAR 的基因组索引。

```bash
STAR --runMode genomeGenerate \
     --genomeDir ~/STAR_data/index \
     --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.107.gtf \
     --runThreadN 4 \
     --genomeChrBinNbits 10 \
     --genomeSAindexNbases 10
```

**参数说明**：

- `--genomeDir`：存放生成的索引的目录。
- `--genomeFastaFiles`：基因组FASTA 文件路径。
- `--sjdbGTFfile`：注释文件 GTF 路径。
- `--runThreadN`：并行线程数（根据CPU核心数调整）。
- `--genomeChrBinNbits`：默认值为 18，将其降低到 10 可以减少内存使用。



**indexing后会得到哪些文件？**

生成的索引文件包括：

- **`Genome`**（二进制编码的基因组序列）
- **`SA`**（后缀数组）
- **`SAindex`**（FM-Index）
- **`sjdbInfo.txt`**（剪切位点数据库）

这些索引文件会在比对过程中被 **STAR** 加载，以加速 reads 的映射。

![Image](https://github.com/user-attachments/assets/16ef8827-80f6-44ee-885e-e3b8ef7e763d)

## 4 RNA-Seq 数据比对

### （1）获取酵母 RNA-Seq 原始数据

我们先获取一组来自酵母 RNA-Seq 的双端配对（paired-end）数据，用于 STAR 测试。

> 数据集介绍：
>
> 数据来源：[https://bioinfogp.cnb.csic.es/files/samples/rnaseq/](https://bioinfogp.cnb.csic.es/files/samples/rnaseq/)
>
> 其提供了一个名为“RNA-Seq_Sample_Files.zip”的数据集，大小为1.3 GB，包含以下内容：
>
> - **测序数据**：多组 Illumina 双端测序的 FASTQ 文件，每个条件下有三个生物学重复（replicate）。
> - **基因组序列**：包含酿酒酵母（Saccharomyces cerevisiae）基因组序列的 FASTA 文件。即为前面下载的酵母基因组序列文件。
> - **基因注释**：包含酿酒酵母基因坐标的 GTF 文件，即为前面下载的酵母基因注释GTF文件。以及一个包含基因符号和描述的制表符分隔的文本文件（ann.txt）。
>
> 这些文件可用于执行一个简单的 RNA-Seq 流程，适合教育和研究目的。
>
> ![Image](https://github.com/user-attachments/assets/03794abd-3094-4ff0-a6eb-84fb4634974c)

**1 创建存放目录**

```bash
mkdir -p ~/STAR_data/test_rnaseq
cd ~/STAR_data/test_rnaseq
```

**2 准备数据**

我们就取3组replicate的双端配对序列作为测试，代表其中的3种实验条件（pH和CO2）：

①1M_SRR9336468_1.fastq.gz（Read 1）、1M_SRR9336468_2.fastq.gz（Read 2）
②1M_SRR9336471_1.fastq.gz（Read 1）、1M_SRR9336471_2.fastq.gz（Read 2）
③1M_SRR9336474_1.fastq.gz（Read 1）、1M_SRR9336474_2.fastq.gz（Read 2）

注：双端测序（Paired-End Sequencing）下，每个片段的两端会分别生成两个配对的读数（reads），存储在两个独立的文件中，通常命名为 **Read 1** 和 **Read 2 ** 或 **File 1** 和 **File 2**。分别代表前向序列（Forward Reads）和反向序列（Reverse Reads）。“SRR9336468”是编号，是一个 SRA 数据库中测序运行的唯一标识符。

> **SRA**（Sequence Read Archive） 是 NCBI 维护的公共高通量测序数据存储库，存储了 Illumina、PacBio、Nanopore 等平台的 RNA-seq 和 DNA-seq 数据。地址：[Home - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra)

我们将这6个gz文件复制到test_rnaseq目录下，解压。

```bash
gunzip *.gz
```



### （2）测试运行 STAR

使用 `STAR` 开始对这组RNA-Seq数据比对。

第①组：

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336468_1.fastq 1M_SRR9336468_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336468_ \
     --outSAMtype BAM SortedByCoordinate
```

**参数说明**

- `--genomeDir`: 之前生成的酵母基因组索引目录。
- `--readFilesIn`: 输入的配对端 FASTQ 文件。输入两段Read。
- `--runThreadN`: 并行线程数。
- `--outFileNamePrefix`: 输出文件的前缀。
- `--outSAMtype`: 输出 BAM 文件格式并按坐标排序。

第②、③组的命令行代码和第①组差不多：

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336471_1.fastq 1M_SRR9336471_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336471_ \
     --outSAMtype BAM SortedByCoordinate
```

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336474_1.fastq 1M_SRR9336474_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336474_ \
     --outSAMtype BAM SortedByCoordinate
```

**输出文件**

运行成功后，会在 `~/STAR_results/` 目录下生成结果文件，如：

- **比对结果（BAM 文件）**: `1M_SRR9336468_Aligned.sortedByCoord.out.bam`
- **比对统计（日志文件）**: `1M_SRR9336468_Log.final.out`



## 5 RNA-Seq 比对结果查看

### （1）查看日志文件的统计数据

执行以下命令查看统计数据：

```bash
cat ~/STAR_results/1M_SRR9336468_Log.final.out
cat ~/STAR_results/1M_SRR9336471_Log.final.out
cat ~/STAR_results/1M_SRR9336474_Log.final.out
```

![Image](https://github.com/user-attachments/assets/48b6a520-85ab-4484-8f60-84d381f2a123)

![Image](https://github.com/user-attachments/assets/dddf9b21-5181-478a-86db-36738669e98d)

![Image](https://github.com/user-attachments/assets/6b7f9394-fa3e-40f7-8745-e979a7479f47)

**以第一个`1M_SRR9336468_Log.final.out`文件为例解读一下：**

输入的reads数（Number of input reads）为 1,000,000
平均读长为（Average input read length）为 300 bp

**唯一比对（Unique Reads）**

> 注：Unique Reads指一个reads可以唯一地比对到参考基因组的一个特定位置。没有任何其他位置能匹配该读数的序列，或者其他位置的匹配分数低于某个阈值。

唯一比对到基因组的读数（Uniquely mapped reads number）：**881,977**；
唯一比对的比例（Uniquely mapped reads）为 **88.20%**。

接着是**剪接**位点的信息（splicing），包括检测到的剪接事件的总数（Number of splices: Total），与GTF文件中注释的剪接位点一致的数量（Number of splices: Annotated (sjdb)），非标准剪接位点（Number of splices: Non-canonical）等。

> 注：剪接位点即基因组中 **外显子 (exon)** 和 **内含子 (intron)** 之间的连接点。**剪接位点数量**表示在测序读数的对齐过程中，STAR检测到的外显子-外显子连接点的数量。
> 关于标准剪接位点 (Canonical Splice Sites)，最常见的剪接模式是 **GT-AG**，即内含子起点为 GT，终点为 AG。
> 非标准剪接位点 (Non-Canonical Splice Sites)常为其他序列模式（如 AT-AC 或 GC-AG）。

**错误率和编辑率**

每个碱基的错配率（Mismatch rate per base）为 **0.25%**，每个碱基的缺失率（Deletion rate per base）为 **0.02%**，每个碱基的插入率（Insertion rate per base）为 **0.01%**。

**多位点比对（Multi-Mapping Reads）**

> 注：Multi-Mapped Reads指一个读数可以比对到参考基因组的多个位置，并且这些位置的匹配分数相近，因此无法明确归属于某一个特定位置。通常来源于参考基因组中具有重复序列的区域。

读数比对到多个基因组位点的数量（Number of reads mapped to multiple loci）为**61,894**，占比 **6.19%**。

**2,484** 个读数比对到太多位点（Number of reads mapped to too many loci），系统将其丢弃。

**未能比对的读数 (Unmapped Reads)**

其中因为太多错配（too many mismatches）而未比对上的读数为 **0**。因为读序太短（too short）而未比对上的读数为 **53,634** ，占比为 **5.36%**。



### （2）查看BAM文件

#### ①使用Linux的Samtools工具 

我们可以使用**Samtools**（最常用的 BAM 文件操作工具，支持快速查看比对统计信息、覆盖深度、特定区域比对等）处理BAM文件。

安装命令：

```bash
sudo apt update
sudo apt install samtools
```

Samtools也可以安装在conda环境里。

**常用命令**：

1. 查看 BAM 文件内容：

   ```bash
   samtools view 1M_SRR9336468_Aligned.sortedByCoord.out.bam | head
   ```

   ![Image](https://github.com/user-attachments/assets/5760a86f-0cf1-4de9-9a00-b9def2095dd2)

2. 查看比对统计：

   ```bash
   samtools flagstat 1M_SRR9336468_Aligned.sortedByCoord.out.bam
   ```

   ![Image](https://github.com/user-attachments/assets/3498b1b7-0b95-4c90-bdf8-ec6339d6da11)

3. 生成覆盖度统计：

   ```bash
   samtools depth 1M_SRR9336468_Aligned.sortedByCoord.out.bam > coverage.txt
   ```

4. 生成BAM文件对应的`.bai`索引文件
   （bai索引的作用是加速访问 BAM 文件中的特定区域，提高数据处理效率。比如常见的可视化基因序列和比对结果的软件`IGV`在载入比对结果时，也需要载入对应的索引文件）

    ```bash
   cd ~/STAR_results/
   samtools index 1M_SRR9336468_Aligned.sortedByCoord.out.bam
    ```



此外，我们可以用IGV软件进行可视化BAM文件，详情：[https://igv.org/doc/desktop/#DownloadPage/](https://igv.org/doc/desktop/#DownloadPage/)



#### ②使用Python库 pysam

**安装**

`pysam`不能跑在Windows上，我们需用conda建立Linux下的Python环境，并安装pysam库。注意，现阶段0.22.1版本的pysam尚不支持Python3.13。

下图展示了在Linux系统的conda环境中安装`pysam`（在WSL的Linux下配置conda环境并使用VS Code详见预备章节1）。

![Image](https://github.com/user-attachments/assets/7ca3fa81-de9a-40dc-b36d-1fc22f55f7b5)

**统计查看**

查看文件内容

```python
import pysam

# 定义 BAM 文件路径
bam_file_path = "../STAR_results/1M_SRR9336468_Aligned.sortedByCoord.out.bam"

# 打开 BAM 文件
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    print(f"读取 BAM 文件：{bam_file_path}")
    print(f"参考基因组名称：{bam_file.header['SQ']}\n")
    
    # 遍历 BAM 文件中的比对
    for i, read in enumerate(bam_file.fetch()):
        print(f"Read {i + 1}:")
        print(f"  QNAME (Query Name): {read.query_name}")
        print(f"  FLAG: {read.flag}")
        print(f"  Reference Name: {read.reference_name}")
        print(f"  Start Position: {read.reference_start}")
        print(f"  Mapping Quality: {read.mapping_quality}")
        print(f"  CIGAR String: {read.cigarstring}")
        print(f"  Sequence: {read.query_sequence}")
        print(f"  Base Qualities: {read.query_qualities}")
        
        # 提取附加标签
        tags = dict(read.get_tags())  # 获取所有标签
        nh = tags.get("NH")  # 比对次数
        hi = tags.get("HI")  # 比对编号
        as_score = tags.get("AS")  # 比对得分
        nm = tags.get("NM")  # 错配数

        # 输出结果
        print(f"  NH (Number of Hits): {nh}")
        print(f"  HI (Hit Index): {hi}")
        print(f"  AS (Alignment Score): {as_score}")
        print(f"  NM (Number of Mismatches): {nm}")
        print()
        
        # 限制输出前 5 条 read
        if i >= 4:
            break
```

![Image](https://github.com/user-attachments/assets/d9516f9c-c9a8-44b4-b3fd-745472b3a74f)

**查看比对统计**

```python
# 初始化统计变量
total_reads = 0
primary_reads = 0 
secondary_reads = 0
supplementary_reads = 0
duplicates = 0
mapped_reads = 0

# 打开 BAM 文件
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    for read in bam_file.fetch(until_eof=True):
        total_reads += 1
        if read.is_secondary:
            secondary_reads += 1
        elif read.is_supplementary:
            supplementary_reads += 1
        else:
            primary_reads += 1

        if read.is_duplicate:
            duplicates += 1

        if read.is_mapped:
            mapped_reads += 1

# 输出结果
print(f"{total_reads} in total (QC-passed reads + QC-failed reads)")
print(f"{primary_reads} primary")
print(f"{secondary_reads} secondary")
print(f"{supplementary_reads} supplementary")
print(f"{duplicates} duplicates")
print(f"{mapped_reads} mapped ({mapped_reads / total_reads * 100:.2f}%)")
```

![Image](https://github.com/user-attachments/assets/f8841d3c-65f8-4e46-8727-562a7f65acf3)

