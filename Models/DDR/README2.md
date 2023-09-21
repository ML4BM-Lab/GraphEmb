#IMPORTANT INFORMATION
## KRON-RLS METHOD (KERNELS)
#https://www.cin.ufpe.br/~acan/kronrlsmkl/

## SUPPLEMENTARY TABLES 
#https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/7/10.1093_bioinformatics_btx731/7/btx731_supplementary_material_ddr.pdf?Expires=1646933866&Signature=c2nsEncs~VbfHA~A7TG~TLBDFFK0mXNkCp5RSp9xqWAtR9VD9Ml4DIq6194O4sNEJF3NWmGAIUeFLcGNUBWaia477WwFiSNnRxntZhf4FiFZz8DiOwyNGbRQSDnGkQppSAxMJLmgpVQhf6bquAiEnONCPOqZ9y3iFcC17GZCHJSm5yjs~uR7Db6FMQVvnZXWzVmb5auVkbG4j5Hh-yDYD2xo~0s0numuVa-ZHgOVsUEz3zKD7O~J9kSywM~msGaSCYBo3~p1~c~VwnS3EEyXLSXlcucGt9V3CQ4ZgUC4ROYRZzz1g6fmgZjzRioT50grwGJgxRMyO6ZlfrSvkaSq3Q__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA


#-------------------------------ADJMATRICES -----------------------------#

#---1---# The gene expression similarities of drugs and of target proteins.

# Data source: Preprocessed CMap files data from (Isik, et al., 2015),
# where we considered only MCF7 cell line instances following the
# method presented in (Hizukuri, et al., 2015).

# Descriptor derivation: Matrix of expression profiles (comprising
# compounds in rows and target proteins in columns), as it explained in
# (Hizukuri, et al., 2015).
# Similarity calculation: The expression similarities of compounds and of
# target proteins, respectively, are calculated by using Pearson’s
# correlation coefficients on the row and column profiles of the expression
# matrix, respectively

# STATUS: WORKING ON IT ...

#---2---# Disease based similarity

# Data source: Drug-disease and target protein-disease associations are
# obtained from KEGG Disease (Kanehisa, et al., 2017).

# Descriptor derivation: Profile of drug–disease pairs and target protein–
# disease pairs that are known to be associated, each drug (or target
# protein) is described by a binary profile represents the presence or
# absence of disease name.
# Similarity calculation: Taminoto coefficient (TC).

# STATUS : WORKING ON IT.. 

#---3---# Pathway based similarity.

# Data source: Drug-pathway and target protein-pathway associations are
# obtained from KEGG Pathways.

# Descriptor derivation: Profile of drug–pathway pairs and target proteinpathway 
# pairs that are known to be associated, each drug (or target
# protein) is described by a binary profile represents the presence or
# absence of pathway name.
# Similarity calculation: TC.

# STATUS : WORKING ON IT.. 

#---4---# Gaussian interaction profile similarity based on the topology of DTI network.

# Data source: Known interactions between drugs and target proteins are
# obtained from DrugBank (Law, et al., 2014).

# It is calculated as in (van Laarhoven, et al., 2011).

# STATUS: WORKING ON IT ...

#---5---# Chemical structures-based molecular fingerprints similarity.

# Data source: Chemical structures are obtained from DrugBank

# Descriptor derivation: Fingerprints (i.e., CDK_Standard, CDK_Graph,
# CDK_Extended, CDK_Hybridization, KR, MACCS, PubChem,
# SIMCOMP, EC4, FC4, EC6, FC6, Lambda, Marginalized,
# MinMaxTanimoto, Tanimoto, and Spectrum) are generated using
# Kebabs (Palme, et al., 2015), Rchemcp (Klambauer, et al., 2015), Rcpi
# (Cao, et al., 2015), CDK (Steinbeck, et al., 2003), and SIMCOMP
# (Hattori, et al., 2010) tools.
# Similarity calculation: TC

# STATUS: WORKING ON IT ...

#---6---# Drug interactions based similarity

# Data source: DrugBank.

# Descriptor derivation: Profile of drug interactions; each drug is
# describing by a binary vector specifying the presence of absence of each
# interacting drug.
# Similarity calculation: TC

# STATUS: WORKING ON IT ...

#---7---# Drug side-effect based similarity

# Data source: SIDER2 (Kuhn, et al., 2016).

# Descriptor derivation: Profile of drug side-effects associations; each
# drug is describing by a binary vector specifying the presence of absence
# of each side effect keyword.
# Similarity calculation: TC

# STATUS: WORKING ON IT

#---8---# Drug ATC-code based similarity

# Data source: DrugBank

# Descriptor derivation: Profile of drug ATC-codes associations; each
# drug is describing by a binary vector specifying the presence of absence
# of each ATC-code.
# Similarity calculation: TC and using similarity-based score from
# (Cheng, et al., 2013).

# STATUS: WORKING ON IT 

#---9---# Protein similarity based on functional annotation using gene ontology (GO)

# Data source: GOA (Barrell, et al., 2009)

# Descriptor derivation: Profile of target protein GO-terms associations,
# for each namespace molecular function (MF), cellular compartment
# (CC) and biological process (BP); each target protein is describing by a
# binary vector specifying the presence of absence of each GO-term.
# Similarity calculation: The semantic similarity is calculated using Rcpi
# tool.

# STATUS: WORKING ON IT ....

#---10---# Protein domain based similarity

# Data source: Pfam (Finn, et al., 2016)

# Descriptor derivation: Profile of target protein-domains associations;
# each target protein is describing by a binary vector specifying the
# presence of absence of each domain.
# Similarity calculation: TC.

# STATUS: WORKING ON IT ...

#---11---# Protein sequence based similarity

# Data source: UniProt, KEGG Genes

# Similarity calculation: It is calculated using a normalized version of the
# Smith-Waterman (SW) algorithm (Smith and Waterman, 1981). We also
# calculated other sequence-based descriptors such as Mismatch and
# Spectrum kernels using Kebabs tool (Palme, et al., 2015).

# STATUS: WORKING ON IT ...

#---12---# Proximity in protein-protein interactions (PPI) network.

# Data source: HIPPIE (Alanis-Lobato, et al., 2017)
# It is calculated as in (Perlman, et al., 2011).

# STATUS: WORKING ON IT ...