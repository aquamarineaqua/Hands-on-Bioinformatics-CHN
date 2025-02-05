## 第一章 RNA-Seq 数据处理与分析：比对、质量控制与定量

**RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification**

本章节讲述 RNA-Seq 数据处理与分析，包括比对、质量控制与定量。本章内容取材自 STAT 115 2021 的 Homework1。

原课程采用人类基因组作为研究对象，并在服务器上跑分析流程，但是普通计算机无法支撑人类基因组的运算。基于教学的目的，为方便普通笔记本也能跑`STAR`等分析，我选择了轻量化的酵母基因组（*Saccharomyces cerevisiae*）作为研究对象。

在实战之前，建议先熟悉原课程中的理论部分。

---

目录：

- **Chapter 1 RNA-Seq 数据处理与分析：比对、质量控制与定量**
  （RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification）
  - **Section 1-1 用 STAR 进行 RNA-seq 测序比对**
  （STAR alignment）
