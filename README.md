# LiBis: An ultrasensitive alignment method for low input bisulfite sequencing

![](https://github.com/Dangertrip/LiBis/blob/master/Figure1.jpg?raw=true)

Methylation information of cell-free DNA in cerebrospinal fluid (CSF), plasma and other body fluids has been utilized to diagnose early-stage diseases and predict therapy response. However, in common clinical settings only a very low amount of cell-free DNA can be purified, usually in the range of a few dozen nanograms. Even worse is that the cfDNA is fragmented and peaked around 120 bases. Whole genome bisulfite sequencing (WGBS), the gold standard to measure DNA methylation, on such a low amount of fragmented DNA molecules introduces a critical data analysis challenge, which is the low mapping ratio. This, in turn, generated low sequencing depth of each CpG and low coverage of genome-wide CpGs sites. The problem of insufficient informative CpGs became the bottleneck of the clinical application of cell-free DNA WGBS assays. Hence, we developed LiBis, a novel method for low input bisulfite sequencing data alignment. By dynamically clipping initially unmapped reads and remapping clipped fragments, we conservatively rescued those reads and uniquely aligned them to the genome. With much improved mapping ratio, LiBis as an integrative toolkit increases the number of informative CpGs and the precision of methylation status of each CpG. High sensitivity and cost efficiency achieved by LiBis for low input samples allow discoveries of genetic and epigenetic features for downstream analysis and biomarker identification from liquid biopsy.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Install Python3, conda and pip
Python3 and conda can be downloaded and installed from https://www.anaconda.com/distribution/
Please make sure than Python version is at least 3.6

Please run following command to make sure python, pip and conda are correctly installed.

```
python --version
pip --version
conda --version
```

#### Install LiBis by Conda (Recommand)

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install libis
```

#### Install LiBis by pip
LiBis integrates FastQC, Bedtools, Trim-galore, moabs and samtools. These packages can be installed independently or by conda.

```
conda install -c bioconda fastqc
conda install -c bioconda bedtools
conda install -c bioconda cutadapt
conda install -c bioconda trim-galore
conda install -c bioconda moabs
conda install -c bioconda samtools=1.1
```    

Check the installation of integrated tools

```
fastqc --version
bedtools --version
trim_galore --version
moabs
samtools --version
```
Install LiBis by pip

```
pip install LiBis
```

#### Check the installation of LiBis

```
LiBis --help
```

The parameters of LiBis will be printed to the screen if install successfully.

#### Docker installation

LiBis also supports the Docker installation.

```
wget https://github.com/Dangertrip/LiBis/archive/master.zip
unzip master.zip
cd LiBis-master/
docker build --tag=libis ./
docker run libis LiBis --help
```

## Running test case

### Run LiBis

```
git clone https://github.com/Dangertrip/LiBis.git
cd Example/

# Run LiBis from begining:
LiBis -n sample1_mate1.fq.gz,sample1_mate2.fq.gz sample2_mate1.fq.gz,sample2_mate2.fq.gz -r /PATH_TO_FASTA_REFERENCE 

```    
We recommand that the original bam files are also mapped by BSMAP.


### Run test case when using Docker(Please put reference and sample raw data under the dictionary you input.)

```
docker run -v /path/to/yourdata:/data/ libis LiBis -f 0 -n name1.fq.gz name2.fq.gz -r hg19.fa
```

## Authors

Yue Yin, Jia Li, Jin Li, Minjung Lee, Sibo Zhao, Linlang Guo, Jianfang Li, Mutian Zhang, Yun Huang, Xiao-Nan Li, Deqiang Sun

## License

This project is licensed under the MIT License.



## Parameters in LiBis

### optional arguments:
### -h, --help            
show this help message and exit
### -f FILE, --file FILE  
Required. Enter a number, 0 means using parameter to set up, 1 means using text file to set up
### -sf SETTINGFILE, --settingfile SETTINGFILE
Required. Setting txt file name. Ignore if -f is 0.
### -n [NAME [NAME ...]], --name [NAME [NAME ...]]
Required. Fastq file name.
### -c CLIP, --clip CLIP  
Clip mode. 0 means close. 1 means open. default: 0
### -l [LABEL [LABEL ...]], --label [LABEL [LABEL ...]]
Required. Labels for samples
### -g GENOME, --genome GENOME
Required. Genome the reference belong to.(Use for plotting) hg18/hg19/mm10/mm9 and so on. Plotting script will not avaliable if leave it blank
### -w WINDOW, --window WINDOW
Window length for clipping mode, default=30
### -s STEP, --step STEP  
Step size for clipping mode.
### -p PROCESS, --process PROCESS
Process using for one pipeline. Normally bsmap will cost 8 cpu number. So total will be 8p.
### -r REF, --ref REF     
Required. Reference
### -qc QUALITYCONTROL, --QualityControl QUALITYCONTROL
Do(1) quality control or not(0)
### -t TRIM, --trim TRIM  
Do(1) trimming or not(0). Don't need to do trimming if you use clip mode.
### -b BINSIZE, --binsize BINSIZE
Plot setting. Set the bin size for averaging methylation ratio among samples, default=1000k
### -ft FILTER, --filter FILTER
filter for clipped reads
### -bam BAM, --bam BAM   
Processed bam file for the first step. If bam files are offered here, the first step of bsmap will be skipped. BAM files can only be generated by BSMAP. Different files should be seperated by ','. If there's no bam file for part of the samples, leave the space blank. For example: a.bam,b.bam,,,,f.bam
### -mcall, --mcall       
Run mcall for mapped bams or not
### -plot, --plot         
Generate the final report or not
### -module, --MOABSmodule
Run LiBis as MOABS module. Please only apply this when integrating LiBis to MOABS.
### -nc, --nocheck        
Skip the checking step for result folders. Using this parameter may rewrite the previous results.
### -bu, --bamu           
Given processed bam file contains unmapped reads. Only use with -bam is open.
### -gz, --gzip           
Temporary files are in gz format.
### -fullmode, --fullmode
Keep all temp files.

## Releasing information

### Version 0.0.1
1. Remove all part fastq/bam/sam. Merge all fragments into one file to speed up the computation.
2. Extension requires the overlap between two fragments.
3. Add -bam option to allow users use their own bam file for the first stage mapping.

### version 0.0.2
1. Decided the name of the software: LiBis, which stands for Low input Bisulfite sequencing 
2. Add -mcall, -plot, -nc, -fullmode
3. normal mode of bsmap now use label as the filename for generated bam.

### version 0.0.3
1. Fixed reads name bug, remove the redundant modification of reads name like "@SRR001666.1". Keep the removal of "/1" or "/2" at the end of the reads.
2. Add left cut length and right cut length to the reads name. Now the reads name contains 4 fields at the end divided by "_":
3. The first number represent the rank of fragment from the reads(For example, if there are two fragments clipped from one reads, the second field will be 0 and 1).  
4. The second one means the mate, which file does the read come from. 
5. The third number means the left cut length
6. The forth number means the right cut length

### version 0.0.9
1. add -module, -mm ,and -pnf.
2. All temp SAM files are replaced by pysam.

### version 0.0.12
1. add -rc. 
2. All temp fastq files are in GZ format now.

### version 0.0.13
1. Remove 'sort' in mapreduce step

### version 0.0.14
1. Add RELEASE to the pypi package

### version 0.0.15
1. Combined multiple diction in mapreduce
2. Reduce the complexity in unmapped reads selection.
3. Changed to order of unmapped reads in unmapped fastq. Now fragments from one unmapped are adjacent and are ranked by file order from 0 to MAX.

### version 0.1.0
1. Add fast mode when BSMAP allows -U to output unsorted reads.
2. MOABS>=1.3.8.5

### version 0.1.1
1. Add -u to LiBis first round mapping.
2. Add reminder for each step.


### version 0.1.2
1. Add 0.05*read_length as mismatch threshold.
2. Now original output bam doesn't contain any unmapped reads.


### version 0.1.3
1. Now use sort in pysam to sort all bam files.
2. Now useless tags will be removed.


### version 0.1.4
1. Now move identifiers in reads' name to tags:
        MT: Mate information. 0: mate1; 1: mate2; 3: single end mapping. 
        FG: Fragment order.
        LC: Clipping point from the reads head.
        RC: CLipping point from the reads tail.
        
### version 0.1.5
1. Now pair end reads in split bam file will be matched to each other. 
2. Fixed bug in clipped bam report and the final report.
