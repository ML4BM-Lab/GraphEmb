###----------------------YAMANISHI-----------------------------###

#To obtain the genes, we can use this page http://rest.kegg.jp/get/{gene}/aaseq
#example -> hsa:10327 -> http://rest.kegg.jp/get/hsa:10327/aaseq

import requests
from lxml import html  
from pubchempy import *

'http://www.genome.jp/dbget-bin/www_bget?drug:'

def getdrug_drugbank(drug):
	xpath_url_2 = '/html/body/main/div/div/div[2]/div[2]/dl[7]/dd[3]/dl/dd[2]/a'
	xpath_url = '/html/body/main/div/div/div[2]/div[2]/dl[7]/dd[2]/dl/dd[2]/a'
	start_url = f'https://go.drugbank.com/drugs/{drug}'
	response = requests.get(start_url)
	tree = html.fromstring(response.text)
	links = tree.xpath(xpath_url)[0].text
	if links:
		return (tree.xpath(xpath_url)[0].text, tree.xpath(xpath_url_2)[0].text)
	else:
		return 'NO_KEGG_ID'


def check_drug_in_kegg(drug):
	r = requests.get(f'http://rest.kegg.jp/list/{drug}').text
	return r

def get_SIMMCOMP_score(drug1, drug2):
	"""
	This function returns the similiraty compound score. 
	The imput is a pair of drugs with the format 'KEGG Drug ID'.
	"""
	score = requests.get(f'http://rest.genome.jp/simcomp2/{drug1}/{drug2}/cutoff=0').text
	return score

def getamino_KEGG(protein):
    r = requests.get(f'http://rest.kegg.jp/get/{protein}/aaseq')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

def getdrug_drugbank(drug):
    r = requests.get(f'https://go.drugbank.com/structures/small_molecule_drugs/{drug}.smiles')
    return r.text

###--------------------------DAVIS ET AL--------------------------------###
#The file to download is http://staff.cs.utu.fi/~aatapa/data/DrugTarget/

#You need to install the api 
#pip install PubChemPy

# from pubchempy import Compound
# import requests

def get_drug_pubchem(drug):
    return Compound.from_cid(drug).isomeric_smiles

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

###----------------------------------DRUGBANK FDA-------------------------------------###
#The drugs can be obtained from the link 
#https://go.drugbank.com/structures/small_molecule_drugs/{DRUG}.smiles
#Example -> Drug: DB00755 -> https://go.drugbank.com/structures/small_molecule_drugs/DB00755.smiles

#The proteins can be obtained from the uniprot link
#https://www.uniprot.org/uniprot/{PROTEIN}.fasta
#Example -> Protein: Q16539 -> https://www.uniprot.org/uniprot/Q16539.fasta


import requests

def getdrug_drugbank(drug):
    r = requests.get(f'https://go.drugbank.com/structures/small_molecule_drugs/{drug}.smiles')
    return r.text

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

###-----------------------------------------------BIOSNAP-------------------------------------------##
#The drugs can be obtained from the link 
#https://go.drugbank.com/structures/small_molecule_drugs/{DRUG}.smiles
#Example -> Drug: DB00755 -> https://go.drugbank.com/structures/small_molecule_drugs/DB00755.smiles

#The proteins can be obtained from the uniprot link
#https://www.uniprot.org/uniprot/{PROTEIN}.fasta
#Example -> Protein: Q16539 -> https://www.uniprot.org/uniprot/Q16539.fasta

import requests

def getdrug_drugbank(drug):
    r = requests.get(f'https://go.drugbank.com/structures/small_molecule_drugs/{drug}.smiles')
    return r.text

def getamino_uniprot(protein):
    r = requests.get(f'https://www.uniprot.org/uniprot/{protein}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return aminoseq

###--------------------------------------------------BINDINGDB--------------------------------------------##
