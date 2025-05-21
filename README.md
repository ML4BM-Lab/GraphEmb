# GraphEmb

This is the official repository for the paper  ***Towards a More Inductive World for Drug Repurposing Approaches***,  
published in _Nature Machine Intelligence_ (February 2025).  📖 [Read the article here!](https://www.nature.com/articles/s42256-025-00987-y)


## 📦 Data & Archival

All data used in the article is openly available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13622942.svg)](https://doi.org/10.5281/zenodo.13622942)

This repository is also archived on Zenodo for long-term preservation and reproducibility: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14068683.svg)](https://doi.org/10.5281/zenodo.14068683)


## 🛠 GUEST Python Package

The **GUEST** Python package used in the study is available here: [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/ML4BM-Lab/GUEST)  
Install it via: ``` pip install graphguest ```

## 📁 How This Repo Is Organized
<p align="center" width="70%">
    <img width="70%" src="https://raw.githubusercontent.com/ubioinformat/GraphEmb/main/imgs/folder_structure.png">
</p>

This repository is divided into 4 blocks:
- **DB**: The databases folder contain all the codes that have been used to preprocess the evaluated main (DrugBank, BIOSNAP, BindingDB, DAVIS and Yamanishi) and complimentary (CTD, FDA, HPRD and SIDER) datasets. Also, there is a link to a zenodo repository to download the available datasets.
- **Models**: Containing a folder for every evaluated model (see Figure 2). Within each folder, all code, input matrices to the model and results is given.
- **N2V**: Containing all the necessary code (and a link to zenodo to download input data to the model) to reproduce the results of every model used, and the figures used in the paper.
- **RMSD_comp**: Containing all the necessary code (and a link to zenodo to download input data) to compute the RMSD matrix.
- **RMSD_validation**: Containing all the necessary code (and a link to zenodo to download input data) to reproduce the *in-sillico* validation performed for Moltrans and HyperAttentionDTI models. 


------

<p align="center" width="85%">
    <img width="85%" src="https://raw.githubusercontent.com/ubioinformat/GraphEmb/main/imgs/graphical_abstract.png">
</p>

To facilitate the benchmark process of different Drug Target Interaction prediction models, we assembled tools and code along this work and presented them as both a comprehensive github repo and a python package, as we believe it can address many limitations of current design and benchmarking processes within the development of *in-sillico* drug repurposing approaches (see Figure 1-G).


The rapid growth of machine learning within the graph embedding and drug repurposing areas has motivated a quick development of methodologies, promoting the fast deprecation of older approaches. This work revealed that some models dependencies have conflicts with recent python packages, or provide uncompleted and unmaintained code. To be able to run all the methods, we *dockerized* them allowing an easy execution on any machine. 

Similarly, building the required complementary matrices for every model is a highly demanding task. Most of the databases used at the moment of the method development are now updated, and this data needs to be retrieved in different ways (e.g., xml, tsv, APIs) with different identifiers (e.g., DrugBank, SIDER, PubChem, Uniprot). Moreover, the code necessary to access this information is often not included within the repositories of the models. For this reason, we provide a GitHub repository that includes the code to generate all the necessary matrices, allowing reproducibility of the drug repurposing evaluation.

Finally, the necessity of graph embedding tools for benchmarking drug repurposing methodologies became evident during the study. For this reason, we developed a Python package (named *GUEST*, available at [https://github.com/ubioinformat/GUEST](https://github.com/ML4BM-Lab/GUEST)) to perform relevant tasks such as splitting the data, testing different criteria or handling negative subsampling. This package can ease the comparison of graph embedding approaches not only for DTIs, but also for other types of graph data. 

