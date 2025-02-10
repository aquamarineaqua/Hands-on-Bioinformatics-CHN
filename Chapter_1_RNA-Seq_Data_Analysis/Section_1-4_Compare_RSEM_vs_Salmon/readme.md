**Chapter 1** RNA-Seq 数据处理与分析：比对、质量控制与定量（RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification）

**Section 1-4** RSEM与Salmon结果对比（Compare RSEM vs Salmon）

本节参考 *STAT 115 2021 Homework 1 Problem 6*

---

我们完成了 RSEM 和 Salmon 这两种不同的 RNA-Seq 定量计算算法。并且在上一节中，我们对比了两种定量计算速度上的差异，Salmon（3秒）远远地快于RSEM（32秒）。

本节我们对比两种定量计算方法在结果上的差异。我们主要通过对比以下几个指标：effective length (有效长度), normalized read counts (标准化读数), TPM 和 FPKM，来分析两者定量结果的差异。

请直接打开 `Compare_RSEM_vs_Salmon(Python).ipynb` 或 `Compare_RSEM_vs_Salmon(R).ipynb`。
