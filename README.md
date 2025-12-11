# **metilene<sup>3**: Identifying DMRs Across Multiple Conditions with Auto-Classification

metilene<sup>3</sup> is a computational tool to identify Differentially Methylated Regions (DMRs) across multiple groups (supervised) or samples (unsupervised). With the identified DMRs, metilene<sup>3</sup> enables inference of epigenetic relationships by constructing a Differentially Methylated Tree (DMTree), which can also be used for sample clustering.

Please see the [metilene<sup>3</sup>-doc](https://zzhu1372.github.io/metilene3-doc) for more details.

![alt text](https://zzhu1372.github.io/metilene3-doc/fig/framework.png "framework")

## Installation
You can download and install metilene<sup>3</sup> on Linux, WSL and macOS from this GitHub repo:
```
git clone https://github.com/zzhu1372/metilene3.git
cd ./metilene3
make
```
Dependencies can be installed with conda:
```
conda create -n metilene3 -c bioconda -c conda-forge python==3.10.0 pandas scikit-learn seaborn biopython gseapy r-base bioconductor-ChIPseeker bioconductor-org.Hs.eg.db bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
conda activate metilene3
```
Please check [here](https://zzhu1372.github.io/metilene3-doc/docs/guide/installation.html) for more details.

## Quick Start
After installation, you can test metilene<sup>3</sup> with the included test dataset ``demo_input.tsv``:
```
python ./metilene3.py -test True -o demo_output
```
You should get these files in the ``demo_output`` folder:
```
DMRs-unsupervised.tsv
DMRs.tsv
clusters.tsv
group-ID.tsv
report.html
```
It should take less than one minute on your laptop.

If success, you can run metilene<sup>3</sup> with your methylation dataset (see [format](https://zzhu1372.github.io/metilene3-doc/docs/guide/quick-example.html#cpg-methylation-matrix-methylationmat)):
```
python ./metilene3.py -i your_methylation.tsv -o your_output
```
Check the [full tutorial](https://zzhu1372.github.io/metilene3-doc/docs/guide) to customize your command. 
