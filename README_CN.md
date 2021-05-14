# rna-seq-pipeline

#### 介绍
转录组分析流程

版本：v0.20(beta)

#### 软件架构
本软件由多个子模块构成：质控、比对、获得基因表达矩阵、基因差异表达分析以及功能富集分析。

目前仅仅支持人类和小鼠数据进行分析


#### 安装教程

    git clone https://gitee.com/chinese_pla_general_hospital/rna-seq-pipeline.git
    Rscript rna-seq-pipeline/install_packages.R

#### 使用说明

本程序不提供参考基因组以及其他数据库文件，因此，请在运行程序之前，确保程序目录下的db数据库文件夹，以及文件夹内的参考基因组、基因组索引文件以及gtf格式的注释文件都存在。

> 在使用本程序之前，请注意安装依赖包：
> 以Ubuntu为例，建议安装`libcurl4-gnutls-dev` `libxml2-dev` `libssl-dev`这几个包，否则可能会导致R语言包安装失败。

    python3 rna-seq-pipeline/wrapper.py \
    -l sample_list.csv \
    -p 2 \
    -o output \
    -i raw_data \
    -s hsa \
    -c compare.csv \
    -j qc,align[Optional. Default is 0, which means won't skip any step.]
    
参考基因组来自ensembl数据库，目前支持仅101版本的注释文件

人类： [人类基因组文件](ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) [基因注释文件](ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz)

小鼠： [小鼠基因组文件](ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz) [基因组注释文件](ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.chr.gtf.gz)

猪：[猪基因组文件](ftp://ftp.ensembl.org/pub/release-101/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz) [基因组注释文件](ftp://ftp.ensembl.org/pub/release-101/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.101.chr.gtf.gz)

下载好基因组和注释文件后请用`gunzip`命令解压到`$local_path/db/Genome/GRC*38/`路径下。对于参考基因组来说，还需要用hisat2工具生成索引文件后才能使用。
我们以人类基因组为例。`$local_path/db/Genome/GRCh38/`路径下应该有两个文件,分别是`Homo_sapiens.GRCh38.dna.primary_assembly.fa`和`Homo_sapiens.GRCh38.101.chr.gtf`
。我们接下来用以下命令生成索引文件：
    
    $local_path/bin/hisat2-2.1.0/hisat2-build \ 
        -p 8 \
        $local_path/db/Genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        $local_path/db/Genome/GRCh38/Homo_sapiens.GRCh38
        
上面的命令中`-p`参数是使用的线程数，接下来是输入的参考基因组文件路径，最后是输出文件的前缀。注意`$local_path`是本程序所在的路径。
运行好上面的命令后，`$local_path/db/Genome/GRCh38/`下应该会多出8个以`.ht2`结尾的文件。参考基因组索引构建完成。

运行之前，我们还需要两个配置文件，分别是"样本信息文件"和"组间比较文件"，这两个文件都是csv格式。

"样本信息文件"格式：

|SampleID|Rep|Group|
|----|----|----|
|Sample1|Rep1|Ctrl|
|Sample2|Rep2|Ctrl|
|Sample3|Rep1|Treat|
|Sample4|Rep2|Treat|
|...|...|...|

第一列是样本名称，第二列是这个样本所对应分组的编号，第三列是样本所在的分组名。注意文件名与原始fastq数据名称要保持一致，比如：
我们有一个样本的两个fastq文件`Sample1_1.fastq`和`Sample1_2.fastq`,"样本信息文件"第一列的样本名就是"Sample1",其中"_1.fastq"和"_2.fastq"
是测序数据的标准后缀，不用把它们写入文件中。

"组间比较文件"格式：

|treat|control|
|----|----|
|Treat|Ctrl|
|...|...|

第一列为实验组组名，第二列为对照组组名。这个文件中的组名必须与样本信息文件内的组名保持一致。



到这里为止，运行前的准备工作全部完成。