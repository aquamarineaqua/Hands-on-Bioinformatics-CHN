## 目录
  - [1 安装相关R包](#1-安装相关r包)
  - [2 R中的数据类型和数据结构](#2-r中的数据类型和数据结构)
  - [3 R向量（Vectors）](#3-r向量vectors)
  - [4 R矩阵（Matrices）](#4-r矩阵matrices)
  - [5 R列表（List）](#5-r列表list)
  - [6 R数据框（Data Frames）](#6-r数据框data-frames)
  - [7 使用`dplyr`进行数据操作](#7-使用dplyr进行数据操作)
    - [7.1 filter (按条件筛选)](#71-filter-按条件筛选)
    - [7.2 select (选择特定列)](#72-select-选择特定列)
    - [7.3 mutate (添加新列)](#73-mutate-添加新列)
    - [7.4 arrange (排序)](#74-arrange-排序)
    - [7.5 summarize (聚合)](#75-summarize-聚合)
    - [7.6 group_by (分组)](#76-group_by-分组)
  - [8 关于数据处理的工作流](#8-关于数据处理的工作流)
  - [9 The Pipe Operator `%>%` (管道操作符)](#9-the-pipe-operator--管道操作符)
    - [与 group_by 一起使用 Pipe Operator](#与-group_by-一起使用-pipe-operator)
    - [练习解答](#练习解答)
  - [10 使用 `ggplot2` 绘图](#10-使用-ggplot2-绘图)
    - [例子①](#例子①)
    - [例子②](#例子②)

本课程参考 Harvard Stat 115/215 Lab1 R Basics 相关内容，作为快速上手Bioinformatics的前置知识，内容简短。

本内容也可直接看对应的ipynb记事本文件。

## 1 安装相关R包

我们在之前的文章中已经从Anaconda创建了R语言环境，且安装了相应的包，并可以用Jupyter Notebook环境写R了。

本节需要用到3个R包：`ggplot2`、`dplyr`、`nycflights13`。

其中ggplot2和dplyr已经包含在r-essentials中了，我们只需要安装nycflights13即可。

```bash
conda activate 环境名
conda install -c conda-forge r-nycflights13
```

注：nycflights13包提供了一个经典的数据集集合，专门用于学习和教学数据分析。它记录了 2013 年美国纽约三大机场（JFK、LGA、EWR）所有航班的数据。可与 dplyr 和 ggplot2 一起用于教学数据清洗和探索性数据分析。

---

## 2 R中的数据类型和数据结构

首先了解R的基本数据类型和数据结构。

R有6种基本数据类型：

- 字符型（character）`"A"`
- 数值型（numeric）`1`
- 整型（integer）`1L`
- 逻辑型（logical）`TRUE, T`
- 复数型（complex）`1+1i`

R有多种数据结构，包括：

- 原子向量（atomic vector）
- 列表（list）
- 矩阵（matrix）
- 数据框（data frame）
- 因子（factors）

> 注：
> 
> ①关于**atomic**（原子性）<br>
> 意为向量或集合中所有元素都是同一种数据类型。
> 
> ②关于**因子**（Factors）数据结构<br>
> **因子**（factors）是R中的一种特殊数据结构，主要用于表示**分类数据**（categorical data）。因子将数据的不同类别存储为**离散的水平（levels）**，并使用整数值来表示这些水平，从而优化内存使用并提高处理效率。

## 3 R向量（Vectors）

R天生支持向量化操作（这与Python中的列表不同）。

```R
# 创建一个向量
# 注意R中对象赋值的符号为：<-
x <- c(1, 2, 3, 4, 5) 
print(x[2]) # 按索引访问元素（索引从1开始）

# 对向量中的每个元素进行操作
x^2
sqrt(x)

# 向量也可以是逻辑型
print(x[x < 3]) # 使用逻辑向量进行索引

# 初始化一个空向量
vector("numeric", 5)
```

![Image](https://github.com/user-attachments/assets/1d840d74-7608-4aba-99d5-b53d88bcad11)

## 4 R矩阵（Matrices）

矩阵是带有维度属性的R向量，它也是原子性的（atomic）。

```R
# 创建一个2x4矩阵
y <- matrix(1:8, nrow = 2, ncol = 4, byrow = FALSE) 
y # 查看矩阵
str(y) # 查看矩阵的结构

# 访问矩阵元素
y[1, 2] # 访问第一行第二列的元素
y[, 2]  # 访问第二列的所有元素
dim(y) # 获取矩阵的维度(size of y)
y %*% t(y) # 矩阵乘法
```

![Image](https://github.com/user-attachments/assets/e962f903-d91e-4667-8382-af9f30cf5cc3)

```R
dim(y) <- NULL  # 将y的维度信息清除，矩阵就变成了展平的一维vector
y
```

![Image](https://github.com/user-attachments/assets/51fd12f6-bf1d-41a3-b0f8-ed29d1cff766)

## 5 R列表（List）

列表是可包含不同类型对象的通用向量。

```R
# 创建一个包含向量、矩阵和列表的列表
list_data <- list(c("Jan", "Feb", "Mar"), 
                  matrix(c(3, 9, 5, 1, -2, 8), nrow = 2),
                  list("green", 12.3))

# 为列表中的元素命名
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# 显示列表的结构
str(list_data)

# 显示列表
list_data
```

![Image](https://github.com/user-attachments/assets/2209ce8b-e74c-4eae-9405-9fb85e78eb90)

## 6 R数据框（Data Frames）

数据框是含有多个相同**长度**的R向量的列表。

```R
data(mtcars) # 加载一个著名的数据集，来自nycflights13
str(mtcars) # 查看数据框的结构，它是数值向量的列表
head(mtcars) # 显示前六行
```

![Image](https://github.com/user-attachments/assets/33f6a2e9-6b6c-4f5b-8d5d-315dc9b9e938)

```R
mtcars[1, 1] # 访问第一行第一列的值
head(mtcars[1]) # 返回数据框中第一列（仍是数据框形式）
head(mtcars[[1]]) # 返回数据框中第一列的向量形式
sapply(mtcars, sum) # 对数据框的每列计算总和
```

![Image](https://github.com/user-attachments/assets/5929c0d0-7399-4605-bd7e-924934d88c67)

## 7 使用`dplyr`进行数据操作

- **`dplyr`** 是一个R包，它常被视为R中一种数据处理的“语言”，可以非常简单直观地操作数据。
- 数据科学（以及计算生物学）的工作中，80%是数据清理，20%是数据分析。
- `dplyr` 的官方文档：[https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)

```R
library(dplyr)
library(nycflights13)

data(flights) # 载入数据集
head(flights)
```

![Image](https://github.com/user-attachments/assets/efd08e3d-f789-42a1-bf64-647ba94376f2)

### 7.1 filter (按条件筛选)

- 根据条件选择特定的行。

```R
# 设置HTML表格显示行数
options(repr.matrix.max.rows = 14) # 设置显示的最大行数
options(repr.matrix.max.cols = 8)  # 设置显示的最大列数（可选）

filter(flights, month == 1, day == 1)
```

![Image](https://github.com/user-attachments/assets/28337fec-183f-448b-bd91-7d92eacc27f9)

### 7.2 select (选择特定列)
- 选择特定的列。

```R
select(flights, year, month, day)
```

![Image](https://github.com/user-attachments/assets/ea8e985c-763b-41b3-a8c5-97aaa09d9bea)

### 7.3 mutate (添加新列)
- 添加新列（可从已有列的数据中创建）。


```R
mutate(flights,
  gain = arr_delay - dep_delay,
  speed = distance / air_time * 60
)
```

![Image](https://github.com/user-attachments/assets/d4a433de-dd3f-4497-af3d-8d94791183af)

### 7.4 arrange (排序)
- 根据列对数据框进行排序。

```R
arrange(flights, year, month, day)
```

![Image](https://github.com/user-attachments/assets/7028a64b-8a57-4faf-b97f-ff0ab7a6ee34)

**desc (降序排序)**

- 使用 desc 按降序排序。


```R
select(arrange(flights, desc(dep_delay)), year, month, day, dep_delay)
```

![Image](https://github.com/user-attachments/assets/bc9adeb3-8904-42c8-9980-ab3f5de3af8f)

### 7.5 summarize (聚合)

- 将多个数字聚合为一个数字。


```R
summarise(flights,
  delay = mean(dep_delay, na.rm = TRUE) # 设定为移除缺失值后再进行计算
)
```

![Image](https://github.com/user-attachments/assets/a5552814-e180-4a70-a462-2238a01a6424)

### 7.6 group_by (分组)

- `dplyr` 的强大之处在于其支持的 `分组-应用-聚合` 工作流。
- 例如，按飞机编号(tailnum)分组，计算航班数量、平均距离和平均延误。


```R
by_tailnum <- group_by(flights, tailnum)
delay <- summarise(by_tailnum,
   count = n(),  # 计算数量
   dist = mean(distance, na.rm = TRUE),
   arr_delay = mean(arr_delay, na.rm = TRUE))
delay <- filter(delay, count > 20, dist < 2000)
delay
```

![Image](https://github.com/user-attachments/assets/84cb1611-70de-4b4a-901d-c25c7d3b61c7)

## 8 关于数据处理的工作流

我们如何一次对数据进行以上多个步骤的处理？请看以下示例：

>```R
>a1 <- group_by(flights, year, month, day)
>a2 <- select(a1, arr_delay, dep_delay)
>a3 <- summarise(a2,
>    arr = mean(arr_delay, na.rm = TRUE),
>    dep = mean(dep_delay, na.rm = TRUE))
>a4 <- filter(a3, arr > 30 | dep > 30)
>```

注意，上面创建了一些变量（如 `a1`），它们仅用于下一步操作后就再也没有被使用。这通常不是一个好的编程习惯，因为浪费了变量。

---

再看这个:

>```R
>filter(
>    summarise(
>        select(
>           group_by(flights, year, month, day),
>           arr_delay, dep_delay
>        ),
>        arr = mean(arr_delay, na.rm = TRUE),
>        dep = mean(dep_delay, na.rm = TRUE)
>  ),
>   arr > 30 | dep > 30
>)
>```

这段代码特别难读！

下面介绍R中最常用的Pipe Operator，它可作为一种标准的数据处理工作流。

---

## 9 The Pipe Operator `%>%` (管道操作符)

如果按“从内到外”嵌套的方式书写代码，代码可能很难读。
为了解决这个问题，可以使用 `%>%` 操作符（Pipe Operator），它可以将 `f(x, y)` 转换为 `x %>% f(y)`的形式。这样代码更容易阅读。

>```R
>flights %>%
>    group_by(year, month, day) %>%
>    select(arr_delay, dep_delay) %>%
>    summarise(
>        arr = mean(arr_delay, na.rm = TRUE),
>        dep = mean(dep_delay, na.rm = TRUE)
>    ) %>%
>  filter(arr > 30 | dep > 30)
>```

Operator我们知道常译作算子，按顺序罗列算子，再根据算子的顺序进行计算，这样的设计非常的舒服和自然。

注：管道操作符 `%>%` 并不是 `dplyr` 发明的，而是来自 `magrittr` 包。它也可以用于其他场景，如：
>```R
>letters %>% length() # 为变量letters计算长度
>```

### 与 group_by 一起使用 Pipe Operator

我们可以与 group_by 一起使用管道操作符

比如以下这段代码：

```R
flights %>%
  filter(origin == 'EWR') %>%
  group_by(dest) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
```

![Image](https://github.com/user-attachments/assets/0a7e6177-d201-4965-b40b-12ae219fe1bb)

**练习：**

我们已经看到的用于数据转换的几个函数：

1) select

2) filter

3) mutate

4) arrange

5) summarise

---

现在需要对 `colon` 数据集按以下步骤进行连续操作（请确保使用 %>% 操作符！）：

1) `study` 列仅包含 1，因此为了简化，删除该列。

2) 我们只需要 etype 为 2 的行（这些是死亡事件，而非复发事件）。

3) 将时间（time）转换为按月计（假设每月正好 30 天）。

4) 按年龄从低到高排序数据框。

5) 按治疗方法（rx）分组，显示一个简单的生存时间（月）的平均值摘要。


---

survival包的`colon`数据集介绍：

>```R
>library(survival)
># 加载数据集
>data(colon)
># 查看数据集的前6行
>head(colon)
>```

`colon` 数据集是 R 的 `survival` 包中包含的一个经典数据集，主要用于研究结肠癌患者的生存分析。该数据集记录了患者的多个临床和试验信息，可用于探索生存时间与各种因素之间的关系。
数据集包含多个变量，包括治疗方法、性别（1 = 男性）、年龄，以及一些试验特定的指标。变量`etype`表示事件类型（死亡或复发）。

---


```R
library(survival)
# 加载数据集
data(colon)
# 查看数据集的前6行
head(colon)
```

![Image](https://github.com/user-attachments/assets/1fba78da-0791-44a5-8452-b60a5f3e355e)

### 练习解答


```R
colon %>%
    select(-study) %>%
    filter(etype == 2) %>%
    mutate(time = round(time / 30)) %>%
    arrange(age) %>%
    group_by(rx) %>%
    summarize(avg_months = mean(time))
```

![Image](https://github.com/user-attachments/assets/8b4a8aef-0958-40b5-b38f-1234689ca8a3)

## 10 使用 `ggplot2` 绘图

我们可以使用`ggplot2`进行绘图。

### 例子①
对 `flights` 数据集进行过滤并创建一个箱线图（Boxplot），显示按月份分类的到达延误情况，并以不同出发地用颜色区分。


```R
library(ggplot2)
flights %>%
    filter(arr_delay <= 360) %>%
    ggplot(aes(x = factor(month), y = arr_delay, color = origin)) +
    geom_boxplot() +
    ggtitle("按月份的延误情况") +
    xlab("月份") +
    ylab("到达延误时间（分钟）")
```

![Image](https://github.com/user-attachments/assets/48a06452-9059-4dcf-9408-770788e3f38b)

ggplot代码详解

**1. 初始化绘图对象：`ggplot(aes(x = factor(month), y = arr_delay, color = origin))`**

- `ggplot()`: 创建一个 ggplot 对象，用于绘图。

- `aes()`

  : 定义图形映射（aesthetics），即图表中的变量和属性的关系：

  - **`x = factor(month)`**: 将 `month` 变量映射到 X 轴，并将其转换为因子（`factor`），以确保月份以分类变量的形式呈现。
  - **`y = arr_delay`**: 将 `arr_delay` 变量映射到 Y 轴，表示到达延误时间。
  - **`color = origin`**: 使用 `origin`（航班出发地）区分数据点的颜色。

**2. 添加几何对象：`geom_boxplot()`**

- `geom_boxplot()`

  : 绘制箱线图，用于显示每个月到达延误时间的分布情况。

  - 每个箱线图表示一个月份的延误数据分布。
  - 通过箱线图可以看到数据的中位数、四分位数范围以及潜在的异常值。

**3. 添加标题和轴标签**

- **`ggtitle("按月份的延误情况")`**: 设置图表的主标题。
- **`xlab("月份")`**: 设置 X 轴标签为“月份”。
- **`ylab("到达延误时间（分钟）")`**: 设置 Y 轴标签为“到达延误时间（分钟）”。

### 例子②

```R
flights %>%
    filter(month == 1, arr_delay < 360, dep_delay < 360) %>%
    ggplot(aes(x = dep_delay, y = arr_delay)) +
    geom_point() +
    ggtitle("Relation between Dep and Arrival Delay") +
    xlab("Departure Delay (min)") +
    ylab("Arrival Delay (min)")
```

![Image](https://github.com/user-attachments/assets/f9e63084-239c-43c4-8cba-5c82f5c61290)
