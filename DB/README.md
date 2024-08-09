Databases
======

You can download the data for running all the models in the zenodo link (X).
Due to size constrains, the datafiles are not uploaded to this GitHub repo.

In most cases, the preprocessed databases are also given.
However, in the ```Code/```  folder you can find the Code needed for preprocessing some of the databases.

## Information about the datasets

In the ```Data/``` folder all the datasets employed in this work can be found in their respective folder:
 - DrugBank. DTIs collected from DrugBank Database Release 5.1.9 [https://go.drugbank.com/](https://go.drugbank.com/). We do not provide the xml files from which we extrated the data because of licensing issues. However, we do provide the DTIs extracted from it, already curated by us. 
 
 - BIOSNAP. Dataset created by Stanford Biomedical Network Dataset Collection. It contains proteins targeted by drugs on the U.S. market from DrugBank release 5.0.0 using MINER. [https://snap.stanford.edu/biodata/](https://snap.stanford.edu/biodata/).

 - BindingDB.  Data employed in this work was the one corresponding to the downloaded from Therapeutics Data Commons (TDC) [https://tdcommons.ai/](https://tdcommons.ai/), however we also provide the data from the original dataset [https://www.bindingdb.org/rwd/bind/index.jsp](https://www.bindingdb.org/rwd/bind/index.jsp). The binarization preprocessing codes are also available.

 - Davis_et_al. Data employed in this work was the one corresponding to the downloaded from Therapeutics Data Commons (TDC) [https://tdcommons.ai/](https://tdcommons.ai/), however we also provide the data from the original dataset. The binarization preprocessing codes are also available.

 - Yamanashi_et_al_GoldStandard. Considered the Gold Standard dataset. It is composed of four subsets of different protein families: enzymes (E), ion channels (IC), G-protein-coupled receptors (GPCR) and nuclear receptors (NR). It was downloaded from [http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/](http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/).
 

Further, in the folder ```cross_side_information_DB``` you can find the datasets used as complementary data to the drug-targer interactions (DTIs) in several models. Inside, you will find the following:

 - CTD. Comparative Toxicogenomics Database for disease-drug and disease-protein associations. Dowloaded from [https://ctdbase.org/](https://ctdbase.org/).

 - FDA. Adverse Event Reporting System (FAERS). The FAERS is a database that contains adverse event reports, medication error reports and product quality complaints resulting in adverse events that were submitted to FDA. Q&A of data: [https://www.fda.gov/drugs/surveillance/questions-and-answers-fdas-adverse-event-reporting-system-faers](https://www.fda.gov/drugs/surveillance/questions-and-answers-fdas-adverse-event-reporting-system-faers).

 - STITCH. Chemical resources from [http://stitch.embl.de/](http://stitch.embl.de/).

 - ChemBL. The file chembl_uniprot_mapping.txt  was used to map between UniprotID and CHEMB identifiers. [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/).
 
 - HPRD. Human Protein Reference Database for human protein-protein interactions. It was downloaded from [http://www.hprd.org/](http://www.hprd.org/).

 - bioMART. It was mainly used to change between Gene Name symbol, UniProtKB, and Gene Name. IT was dowloaded from [https://www.ensembl.org/info/data/biomart/index.html](https://www.ensembl.org/info/data/biomart/index.html).

 - DrugBank. DrugBank Database can be used for extracting other information. In this folder we provide the corresponding KEGG identifiers for DB compounds. 

 - DrugBank_FDA. Dataset of Drugbank compounds approved by FDA.

 - SIDER. Side Effect Resource Database aggregates information from side effects. Dowloaded from [http://sideeffects.embl.de/](http://sideeffects.embl.de/)

 
			