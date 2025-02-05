Windows下部署Linux系统和开发环境

## 目录
  - [1 通过 WSL 安装Ubuntu系统](#1-通过-wsl-安装ubuntu系统)
  - [2 在 WSL 中安装 conda](#2-在-wsl-中安装-conda)
  - [3 Linux下用conda管理环境，并在VS Code中使用环境](#3-linux下用conda管理环境并在vs-code中使用环境)

我们进行生物信息学分析时，几乎所有常用的工具都是安装在Linux系统下，以命令行的形式进行数据处理、分析的。此外，Python的一些用于生物信息学分析的库也是对于Linux支持更好。

所以如果你使用的是Windows系统，则需要想办法让你的计算机也能跑起Linux系统。除了传统的安装双系统和使用虚拟机的方法之外，现在Windows提供了**WSL**（Windows Subsystem for Linux）功能，即Windows子系统，让Windows计算机可同时访问 Windows 和 Linux。

本文讲解如何**用WSL在Windows下部署Ubuntu系统**（最常用的Linux系统），并部署好开发环境。


## 1 通过 WSL 安装Ubuntu系统

确保你的Windows系统是最新版本（11或10的最新版本）。然后打开**命令提示符**（或**PowerShell**），输入：

```powershell
wsl --install
```

该命令会执行一系列操作，包括下载并安装最新的 Ubuntu Linux 发行版（如图）。

有兴趣的话，更多信息可以参考官方教程：①[https://learn.microsoft.com/zh-cn/windows/wsl/install](https://learn.microsoft.com/zh-cn/windows/wsl/install)、②[https://learn.microsoft.com/zh-cn/windows/wsl/setup/environment](https://learn.microsoft.com/zh-cn/windows/wsl/setup/environment)

![Image](https://github.com/user-attachments/assets/793d1101-62c5-4a7b-8579-7285ed009a6d)

提示已安装Ubuntu后，重启Windows，Ubuntu会添加到开始菜单中，我们打开Ubuntu图标。

![Image](https://github.com/user-attachments/assets/c927211c-22e8-42b1-839b-bfaa219494ad)

稍等一会儿后，设置Linux系统的用户名和密码，成功进入Ubuntu系统。

可以看到系统自动安装了最新的Ubuntu24.04版本。这个界面就是Linux系统的命令行界面。

![Image](https://github.com/user-attachments/assets/b58665ee-d198-454b-848d-acb8497c43f7)

第一次进入系统可以使用以下命令更新和升级一下系统内的软件包：

```bash
sudo apt update && sudo apt upgrade
```

---

现在我们可以在Windows中访问Ubuntu的文件。Ubuntu系统的文件地址一般为：`\\wsl.localhost\Ubuntu`，我们也可以在资源管理器中打开它。

![Image](https://github.com/user-attachments/assets/68bf56db-d770-459d-b667-d36895a817d3)

此外，在Ubuntu的命令行界面中，可使用Windows文件管理器打开当前Ubuntu目录。命令如下（注意最后是空格+一点）：

```bash
explorer.exe .
```

注：更多文件系统共享可以参考官方文章：[https://learn.microsoft.com/zh-cn/windows/wsl/filesystems#file-storage-and-performance-across-file-systems](https://learn.microsoft.com/zh-cn/windows/wsl/filesystems#file-storage-and-performance-across-file-systems)

---

我们同样可以使用 VS Code 进行编程。

进入目标目录后，输入`code .`，选择允许主机即可进入VS Code。

![Image](https://github.com/user-attachments/assets/9b6db745-ac1d-47e9-b593-51e60f353bec)

进入后，我们可以新建一个`ipynb`记事本文件，然后保存，保存的时候如果找不到Linux的目录，可以打开文件管理器将当前目录的路径复制下来：

![Image](https://github.com/user-attachments/assets/55487394-3408-4954-a4cb-40526f803020)

然后在VS Code的文件保存界面，粘贴该路径即可：

![Image](https://github.com/user-attachments/assets/0be154c7-462f-48f1-b591-2b8ea4d1336f)

更多WSL上的VS Code相关教程：[https://learn.microsoft.com/zh-cn/windows/wsl/tutorials/wsl-vscode](https://learn.microsoft.com/zh-cn/windows/wsl/tutorials/wsl-vscode)

以上可在Windows上的Python环境中进行Coding，下文将讲解在Linux中用`conda`进行环境管理，并让VS Code识别并使用 WSL的Linux 中的conda环境以进行编程等操作。

## 2 在 WSL 中安装 conda

由于Windows上的conda环境和WSL中conda环境不通，并且我们需要在Linux上使用各类生物信息学工具，所以在 WSL 中安装 Linux 版的 Miniconda 或 Anaconda 是必要的。这里我们选择安装 `Miniconda`。

有了conda工具，我们就可以在Linux中管理Python环境和程序库（如各种生信分析工具`STAR`,`samtools`等）了。

**1 下载 Miniconda for Linux**
在 Ubuntu 的终端中运行以下命令（安装命令可以在Anaconda的官网查询）：

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

**2 运行安装脚本**

安装文件下载到当前目录后，直接在当前目录运行进行安装。

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

**3 按提示操作并安装Miniconda**

- 阅读并同意许可协议。
- 安装路径建议使用默认路径：`/home/<your-username>/miniconda3`。

**4 将 Conda 添加到 PATH**

如果安装完成后仍然无法运行 `conda` 命令，可以手动将 Conda 的路径添加到环境变量。

使用文件管理器直接打开 `.bashrc` 文件。

![Image](https://github.com/user-attachments/assets/c0659f4b-7538-4743-b27d-40992a5db123)

添加下面这一行：

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
```

![Image](https://github.com/user-attachments/assets/0d4a9fe9-03f8-43ef-8750-0ff197127d71)

保存。然后在终端输入：

```bash
source ~/.bashrc
```

![Image](https://github.com/user-attachments/assets/435fa7b3-c4d5-4c34-9ea6-01594b88b40c)

之后输入 `conda --version`，如果正常显示版本号，即成功添加环境变量。

之后我们关闭并重启一下终端即可。

## 3 Linux下用conda管理环境，并在VS Code中使用环境

现在我们可以在WSL的Linux中用`conda`进行环境管理了，接下去我们就可以在 Windows 下用 VS Code，识别并使用 Linux 的 conda 环境进行编程等操作。

为什么要这样做呢？如前所述，因为Linux的conda环境，既能管理Linux的各类应用程序（如生信分析工具`STAR`,`samtools`等），又方便管理Python库，并且有些库如`pysam`（用于python解析STAR的BAM文件）只支持Linux。

---

例如，我们可以在Linux中运行以下命令来创建一个conda环境，既安装对应的Python、R和他们的常用包，又安装想要的程序。

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

之后我们就可以在Linux的这个新环境 r_bio 下进行各项生信分析作业了。

---

接下来我们让 VS Code 能识别到 WSL的Linux 中的各个 conda 环境。

打开VS Code，进入插件商店，搜索`Remote Development`，进行安装。

![Image](https://github.com/user-attachments/assets/afebf689-61d0-4287-9b44-16a89c3095cf)

安装完后，我们在VS Code左侧的远程资源管理器中可以选择WSL作为目标了。这样我们就可以选择使用Linux的 conda 环境了，VS Code会自动探测到已创建的环境，见下图。（注意我们需要为WSL环境再安装一遍必要的插件，比如Jupyter、Python、R等。直接在VS Code的插件界面安装即可）

![Image](https://github.com/user-attachments/assets/ef2f33c0-3d8b-4d13-89a9-1208a51790a7)
