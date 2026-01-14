# Deconomix-case-study
An example case study for the Deconomix python package. It demonstrates the capabilities of the package in a breast cancer scenario, including data curation, preprocessing, processing and visualizations.

# Prerequisites
While the analysis part in this study is computable on standard personal computer, the preprocessing of the data is not feasible on such a system. Therefore we strongly advise you to run this study on a cpu cluster with around 200GB of RAM, which will also drastically reduce the compute time for the hyperparameter search, where a multitude of models have to be trained. We provided an example `Dockerfile` as well.

# Install instructions
To run this case study on your own system, you can clone this repository to your disk. The source of the python package is included as a git submodule, therefore you should run following commands to run the analysis:

```
git clone https://github.com/Deconomix/Deconomix-case-study.git
cd Deconomix-case-study
git submodule update --init --recursive
```

Then you can create a Python virtual environment or a Docker container and run this command inside, which will then install Deconomix form the git submodule and all other requirements needed for data curation and preprocessing.

```
pip install -r requirements.txt
```
# Dataset
In this study two public datasets will be used that have to be downloaded separately:

 - DISCO Breast Cancer/Healthy Single-Cell Expression Data: Download `disco_breast_v01.h5ad` from https://zenodo.org/records/7396984
 - TCGA Breast Cancer Bulk Expression Data. Download via `tcga_downloader.py` which implements an API call to the GDC database

 Move these files to the `Data` directory and then all scripts should be able to run.