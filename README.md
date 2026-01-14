# Deconomix-case-study
An example case study for the Deconomix python package. It demonstrates the capabilities of the package in a breast cancer scenario, including data curation, preprocessing, processing and visualizations.

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
