# **metilene<sup>3**: Identifying DMRs Across Multiple Conditions with Auto-Classification

metilene<sup>3</sup> is a computational tool to identify Differentially Methylated Regions (DMRs) across multiple groups (supervised) or samples (unsupervised). With the identified DMRs, metilene<sup>3</sup> enables inference of epigenetic relationships by constructing a Differentially Methylated Tree (DMTree), which can also be used for sample clustering.

Please see the [metilene<sup>3</sup>-doc](https://zzhu1372.github.io/metilene3-doc) for more details.

![alt text](https://zzhu1372.github.io/metilene3-doc/fig/framework.png "framework")

## Installation
```bash
git clone https://github.com/zzhu1372/metilene3.git
cd ./metilene3
```

### Option 1: Using Pixi (Recommended for development)
If you use [Pixi](https://pixi.sh), you can set up the environment and run tasks instantly:
```bash
pixi install
pixi run test
```

### Option 2: Using containers (Docker & Apptainer)

#### Docker
Build the Docker image locally from the repository root:
```bash
docker build -t metilene3:latest .
# Run metilene3 using Docker:
docker run --rm -v $(pwd):/data metilene3:latest -i demo_input.tsv -g group_info.tsv -o demo_output
```

#### Apptainer / Singularity
Build the Singularity image from the provided recipe:
```
apptainer build metilene3.sif Singularity
# Run metilene3 using Apptainer:
./metilene3.sif -i demo_input.tsv -g group_info.tsv -o demo_output
```

### Option 3: Using Make & Conda / Mamba
Compile locally using make:
```bash
make
```
Dependencies can be installed with conda:
```
conda create -y -n metilene3 -c bioconda -c conda-forge python==3.10.0 pandas numpy matplotlib pandarallel scikit-learn seaborn biopython gseapy r-base bioconductor-ChIPseeker bioconductor-org.Hs.eg.db bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
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
