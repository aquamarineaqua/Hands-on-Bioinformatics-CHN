
使用Anaconda管理R语言环境，并使用Jupyter Notebook编写R语言

## 目录

  - [Anaconda中创建R环境](#anaconda中创建r环境)
    - [0 官方教程存在的问题](#0-官方教程存在的问题)
    - [1 创建R语言环境](#1-创建r语言环境)
    - [2 安装常用包集合：r-essentials](#2-安装常用包集合r-essentials)
    - [3 用VS Code的Jupyter插件写R](#3-用vs-code的jupyter插件写r)




相信一直使用Python做数据分析、机器学习的同学们会习惯使用 Anaconda管理不同的Python环境，并使用Jupyter Notebook（包括使用VS Code的Jupyter插件）来编写Python代码。

切换到R语言时，我们并不习惯使用RStudio这一R常用的代码环境。第一，官方的R语言+RStudio的组合默认是在全局环境下安装管理R包的，没有像Anaconda那样方便管理不同R语言环境的功能。第二，RStudio进行R语言代码编写类似直接编写py文件并运行，不像Jupyter Notebook一个个代码块分别运行和显示那样直观。

本文教大家使用Anaconda来管理R语言环境，并使用Jupyter Notebook编写R程序。本教程需要先安装Anaconda或Miniconda。本教程Windows和Lunix通用。

## Anaconda中创建R环境
### 0 官方教程存在的问题
使用Anaconda官方的教程（[https://docs.anaconda.com/working-with-conda/packages/using-r-language/](https://docs.anaconda.com/working-with-conda/packages/using-r-language/)）创建R环境（`conda create -n r_env r-essentials r-base`）会有一个问题，创建R的最高版本只有3.6。这是因为官方的channel未包含新的R版本。
可以在终端通过命令`conda search r-base`来看所有的r-base的版本：

![Image](https://github.com/user-attachments/assets/c662bfc8-ec0f-4eb3-b155-405950f5781d)

可以看到最新R语言版本的支持的channel位于`conda-forge`中。而`pkg/r`中对R语言的支持没有更新到新版本。

### 1 创建R语言环境

因此，我们需要通过**conda-forge**这一channel安装最新的R语言版本环境。

打开Anaconda Prompt或命令提示符。

我们可以先创建一个名为“*r_ds*”的环境，可采用R语言版本为4.4.1，从conda-forge安装。

```bash
conda create -n r_ds -c conda-forge r-base=4.4.1
```

**注意**，这里我们可以同时把Python的版本也在这里指定好，这样顺带把Python环境也一起配置好了。如果这里不指定这一步，Python会在下面安装r-essentials的时候附带一起安装，这样会默认安装最新的Python版本。但如果你日后常用一些Lunix下的生物信息学分析工具，还是不推荐自动安装最新版本的Python，我这里指定Python版本为3.10。

代码如下：

```bash
conda create -n r_bio -c conda-forge r-base=4.4.1 python=3.10
```

如图安装成功！

![Image](https://github.com/user-attachments/assets/c4eadaad-9862-4606-b315-b348ce757ac2)

### 2 安装常用包集合：r-essentials

创建完我们想要的R环境后，我们可以为其安装单独的R包，可以直接使用命令`conda install -c conda-forge r-包名称`安装。基本上所有常见的R包在Anaconda的环境下都是以`r-`开头的。

我们也可以选择直接安装R基础包集合（R Essentials bundle），即r-essentials，里面包含了80多个常见的R包扩展，如IRKernel, dplyr, shiny, ggplot2, tidyr, caret, nnet等。

注：anaconda中的R包详情详见：[https://repo.anaconda.com/pkgs/r/](https://repo.anaconda.com/pkgs/r/)。 日后可以用来查找某个R包是否在这里面，如果存在的话就可以直接用上面提到的`conda install -c conda-forge r-包名称`命令来安装了。

我们先激活环境，然后安装：

```bash
conda activate r_ds
conda install -c conda-forge r-essentials
```

r-essentials安装完毕后，我们的Jupyter将会支持R内核（如下图），我们在终端输入`jupyter notebook`即可打开记事本界面。这样我们就可以在终端当前的根目录下进行创建ipynb记事本文件并使用jupyter的环境进行R语言编写了。

![Image](https://github.com/user-attachments/assets/04336901-1e36-4a5b-9f18-8089c1d54be0)

---

以上安装成功后，我们就可以从Jupyter notebook启动R语言终端了。终端输入`jupyter notebook`即可。

![Image](https://github.com/user-attachments/assets/afe63b7a-99a9-417e-a8f9-dd177f646985)

![Image](https://github.com/user-attachments/assets/8c0b4512-e93f-4dba-959f-72e184e4bce0)

### 3 用VS Code的Jupyter插件写R

对于平时常用VS Code的Jupyter插件写Python的同学，一定也想用同样的方法进行R语言的编写。前面我们已经从Anaconda安装好R环境`r_ds`了，下面是VS Code中的步骤。

①打开VS Code，如图安装R语言插件。（注意VS Code的Jupyter插件要先安装好）

![Image](https://github.com/user-attachments/assets/46fdceed-ee95-4a4f-b06d-37450a1acc62)

②创建ipynb文件，进入后单击右上角的环境选择按钮，选择环境。其中选择“Jupyter Kernel”。

![Image](https://github.com/user-attachments/assets/2b3d4da8-d61f-4831-81a4-9e5347f42880)

之后选择你刚才从Anaconda创建的R环境，比如`r_ds`。这样，我们就可以成功地在VS Code的Jupyter Notebook中进行R语言编写了！

![Image](https://github.com/user-attachments/assets/31b8ec68-501e-46f8-8b6e-b9b896a5b309)

---

展示下我使用VS Code的Jupyter Notebook写R的示例，是不是很舒适~

![Image](https://github.com/user-attachments/assets/692a9b24-da79-43b5-b146-613b3111dc7a)
