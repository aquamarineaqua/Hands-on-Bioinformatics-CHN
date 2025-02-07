**Chapter 1** RNA-Seq 数据处理与分析：比对、质量控制与定量（RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification）

**Section 1-3** RNA-Seq 定量（RNA-Seq quantification）

本节参考 *STAT 115 2021 Homework 1 Problem 4,5*

## 目录
  - [0 准备：创建新环境](#0-准备创建新环境)
  - [1 RNA-Seq定量（RNA-seq quantification）](#1-rna-seq定量rna-seq-quantification)
    - [（1）安装RSEM](#1安装rsem)
    - [（2）STAR+RSEM](#2starrsem)
    - [（3）Salmon](#3salmon)
      - [①安装Salmon](#①安装salmon)
      - [②Salmon索引目录生成](#②salmon索引目录生成)
      - [③使用双端测序的FASTQ数据运行 Salmon](#③使用双端测序的fastq数据运行-salmon)
    - [（4）对比 STAR+RSEM 和 Salmon 的速度。](#4对比-starrsem-和-salmon-的速度)



## 0 准备：创建新环境

我们在之前的实战中使用`conda`创建过一个新环境 *star_env*，并在里面安装了STAR、samtools等程序。

为了既方便我们使用STAR等Linux程序，又方便我们使用Python和R，现在起我们创建一个新的环境（如命名为`r_bio`）整合一下所有的功能。我们在其中安装R和Python环境（方便我们用VS Code等写程序），同时安装STAR、samtools等程序。这样我们就不用老是切换环境了。

bash命令集如下，我们选择 `r-base=4.4.1`，`python=3.10`。
该指令已经在“预备章节3”中提及，如果你已经创建了这么一个环境，则跳过本小节即可。

```bash
conda config --add channels defaults
conda config --add channels conda-forge  # 添加channel
conda config --add channels bioconda  # 添加channel

conda create -n r_bio -c conda-forge r-base=4.4.1 python=3.10  # 创建新环境r_bio，指定Python版本和R版本

conda activate r_bio  # 激活进入环境

conda install -c conda-forge r-essentials  # 安装常用R包

conda install numpy pandas matplotlib openpyxl scipy sympy jupyter # 安装常用Python库
conda install biopython
conda install pysam

conda install -c bioconda star  # 安装生物信息学序列比对分析程序STAR
conda install -c bioconda samtools  # 安装生物信息学数据处理程序samtools
```




## 1 RNA-Seq定量（RNA-seq quantification）

转录本定量（Transcript quantification）在 RNA-seq 数据分析中起着重要作用，目前已有大量工具用于在转录本层面的量化表达。**`RSEM`**（*Bo Li et al, BMC Bioinformatics 2011*）软件包可从单端或双端 RNA-Seq 数据中估算 *gene* 和 *isoform* 的表达水平。RSEM它支持3种不同的比对工具：bowtie、bowtie2 或 STAR。

**`Salmon`**（*Rob Patro et al, Nature Methods 2017*）是一种超快速、**无需比对**的表达定量的方法，同时能够校正 GC 偏倚。

下面我们选择前面的 **1M_SRR9336471** 这个比对质量高的样本，对其分别使用**STAR+RSEM** 和 **Salmon** 进行定量计算。

### （1）安装RSEM

```bash
conda activate r_bio
conda install -c bioconda rsem
```

### （2）STAR+RSEM

**①STAR**

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn ~/STAR_data/test_rnaseq/1M_SRR9336471_1.fastq ~/STAR_data/test_rnaseq/1M_SRR9336471_2.fastq \
     --runThreadN 4 \
     --quantMode TranscriptomeSAM \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336471_ \
     --outSAMtype BAM SortedByCoordinate
```

关于参数：`--quantMode TranscriptomeSAM` 

默认情况下，STAR 生成的 BAM 文件是基于 **基因组（genome）** 坐标的比对结果。加上 `--quantMode TranscriptomeSAM`，STAR 还会额外输出一个基于 **转录组（transcriptome）**坐标 的 BAM 文件（`Aligned.toTranscriptome.out.bam`）。这个文件仅包含唯一比对（uniquely mapped reads），并将 reads 映射到参考转录本上，而不是基因组坐标。之后，该文件**为后续定量分析提供输入**，适合如 RSEM、Salmon 或 kallisto 用于转录本层面的定量分析。



**②RSEM索引目录生成**

我们首先需要通过RSEM命令行工具 `rsem-prepare-reference` 生成RSEM索引目录。它将参考基因组和基因注释文件（FASTA 和 GTF 文件）预处理为 RSEM 可用的索引文件，用于后续基因和转录本的表达定量。

创建新文件夹，存放RSEM索引。

```bash
mkdir -p ~/RSEM_data/index
```

生成 RSEM 索引。指定GTF格式的注释文件，基因组的FASTA文件，输出目录。（注意命令行中'rsem_index'是输出文件的前缀名）。

```bash
rsem-prepare-reference --gtf ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.gtf \
~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
~/RSEM_data/index/rsem_index
```



**③RSEM定量计算**

先新建目录以存放结果：

```bash
mkdir -p ~/RSEM_results
```

---

命令行如下：

我们将之前STAR中得到`1M_SRR9336471_Aligned.toTranscriptome.out.bam`文件作为RSEM的输入；`~/RSEM_data/index/rsem_index`作为index。注意，由于我们使用了双端测序（paired-end RNA-Seq），我们需要在命令中明确指定 `--paired-end` 参数，否则 RSEM 默认将输入文件视为单端（single-end）数据。结果存储在`~/RSEM_results/RSEMOut`中。

```bash
rsem-calculate-expression --no-bam-output \
--paired-end \
--time \
--bam \
-p 8 \
~/STAR_results/1M_SRR9336471_Aligned.toTranscriptome.out.bam \
~/RSEM_data/index/rsem_index \
~/RSEM_results/RSEMOut
```

可以看到生成的`RSEMOut.isoforms.results`文件中包含了每一个transcript的length, effective length, expected  count, TPM, FPKM, IsoPct信息。

![Image](https://github.com/user-attachments/assets/67be8167-039d-4dc2-af29-782f92421b0a)

### （3）Salmon

#### ①安装Salmon

```bash
conda activate r_bio
conda install -c bioconda salmon
```

#### ②Salmon索引目录生成

**（1）生成转录本序列文件**

使用Salmon做定量计算，首先需要通过salmon命令行工具 `salmon index` 生成Salmon的**转录本索引**。我们需要准备好FASTA格式的 **转录本序列文件**（Transcriptome Sequence），用于后续分析。

注意，**基因组序列（Genome Sequence）** 和 **注释文件（Annotation File）** 可以结合生成 **转录本序列（Transcriptome Sequence）**，因为注释文件定义了基因组中每个转录本的结构（如外显子的位置、顺序和拼接方式），而基因组序列提供了实际的 DNA 碱基序列。

因此，我们需要先从基因组序列文件（`Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa`）和注释文件（`Saccharomyces_cerevisiae.R64-1-1.107.gtf`）生成转录本序列文件。

我们用 **`gffread`** 工具进行从注释文件+基因组序列到转录本序列的转换（[具体可见官网文档](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)）。安装gffread：

```bash
conda install -c bioconda gffread
```

从GTF 文件和基因组序列生成转录本序列：

```bash
mkdir -p ~/Salmon_data/transcriptome_sequence
gffread ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.gtf \
        -g ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -w ~/Salmon_data/transcriptome_sequence/transcripts.fa
```



**（2）生成Salmon用的转录本索引**

有了转录本序列后，可以生成转录本索引。

创建新文件夹，存放Salmon索引。

```bash
mkdir -p ~/Salmon_data/index
```

生成 Salmon 索引。`-t`指定转录本序列文件，`-i`指定输出的索引文件目录，`-k`为k-mer 的长度（默认 31，可根据读长调整）。

```bash
salmon index -t ~/Salmon_data/transcriptome_sequence/transcripts.fa \
-i ~/Salmon_data/index \
-k 31
```

> 关于 Decoy Sequence 警告
>
> 生成索引时可能会提示：
>
> *The salmon index is being built without any decoy sequences. It is recommended that decoy sequence (either computed auxiliary decoy sequence or the genome of the organism) be provided during indexing.*
>
> 这是因为 Salmon 建议在生成索引时包含 **decoy sequences（诱捕序列）**，通常是参考基因组序列或其他非转录区域的序列。Decoy sequences 的作用是减少虚假比对（mapping ambiguity），特别是在存在非特异性 reads 时，可以提高比对的准确性。如果不添加 decoy sequences，Salmon 的定量可能会受到一些非特异性比对的干扰，尤其是对于存在高重复序列或基因组复杂的物种。
>
> 具体可参考：[https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode)
>
> 这里我们暂先不考虑decoy sequences。



#### ③使用双端测序的FASTQ数据运行 Salmon

```bash
salmon quant -i ~/Salmon_data/index \
             -l A \
             -1 ~/STAR_data/test_rnaseq/1M_SRR9336471_1.fastq \
             -2 ~/STAR_data/test_rnaseq/1M_SRR9336471_2.fastq \
             -p 4 \
             -o ~/Salmon_data/output
```

参数说明：
`-l A`：自动检测读向。
`-1` 和 `-2`：分别指定双端测序的 R1 和 R2 文件（如果是单端测序则用参数`-r`）。
我们继续选择使用 1M_SRR9336471 这一双端测序样本。
`-p 4`：使用的线程数。
`-o`：输出结果目录。

---

运行完成后，在输出目录output下的 `quant.sf` 文件中，可以看到包含了 Length、EffectiveLength、TPM、NumReads等信息。

![Image](https://github.com/user-attachments/assets/7923cf07-bcc1-469e-81fa-69ab846831c5)

### （4）对比 STAR+RSEM 和 Salmon 的速度。

通过查询日志文件：

**Salmon** (`Salmon_data/output/logs/salmon_quant.log`):

```
start_time: [2025-01-28 15:24:59.905]
end_time: [2025-01-28 15:25:02.840]
```

**STAR** (`STAR_results/1M_SRR9336471_Log.final.out`):

                                 Started job on |	Jan 27 14:39:06
                             Started mapping on |	Jan 27 14:39:06
                                    Finished on |	Jan 27 14:39:24

**RSEM** (`RSEM_results/RSEMOut.time`)

```
Aligning reads: 0 s.
Estimating expression levels: 14 s.
Calculating credibility intervals: 0 s.
```

---

- Salmon: 大约3秒

- STAR+RSEM: 18秒+14秒=32秒
