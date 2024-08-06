DDR
====

# https://bitbucket.org/RSO24/ddr/src/master/
# 
# <----YAMANISHI----->

# not yet

## DRUGS <---#

Side Effect /*
#    #AERSBIT ^^ FDA, for those not found -> 0 (Weighted cosine correlation coefficient) ✓✓ (IN MARGARET)

#    #AERSFREQ ^^ FDA, for those not found -> 0 (Weighted cosine correlation coefficient) ✓✓ (IN MARGARET)
*/

Drug Properties /*
#    #LAMBDA ^^ RchemCPP --> WORKING ON IT (sd2gramSpectrum) ✓✓

#    #MARGINALIZED ^^ RchemCPP --> WORKING ON IT (sd2gram or sd2gramSpectrum) ✓✓

#    #MINMAXTANIMOTO ^^ RchemCPP --> WORKING ON IT (sd2gramSpectrum) ✓✓

#    #SPECTRUM ^^ RchemCPP --> WORKING ON IT (sd2gramSpectrum) ✓✓

#    #TANIMOTO ^^ RchemCPP --> WORKING ON IT (sd2gramSpectrum) ✓✓
*/

Side Effect
#    #SIDER ^^ Presence or abscence 1 or 0, not found -> 0 (Weighted cosine correlation coefficient) ✓✓/2 (IN MARGARET)

Drug Properties
#    #SIMCOMP ^^ (Yamanishi et al) --> ✓ (IN MARGARET)

SMILES ARE NOT USED EXPLICITLY

## PROTEINS <--#

#    #GO ** BioMART DB, R package -> csbl.go ✓

#    #MISMATCH (4 levels) ** R Package -> KeBABS ✓✓

#    #SPECTRUM (2 levels) ** R Package -> KeBABS ✓✓

#    #PPI ** BioGRID DB, using paper formula (is in his matlab scripts, 'shortest hop distance') ✓✓

#    #SW ** ALREADY HAVE INFO FOR THAT --> ✓✓ (IN MARGARET)


# ---------------------------------------------------------------------#
# Entity Kernels Information source
# -----------------------------------
# Drugs       AERS-bit - AERS bit Side-effects
#             AERS-freq - AERS freq Side-effects
#             GIP - Gaussian Interaction Profile Network
#             LAMBDA - Lambda-k Kernel Chem. Struct.
#             MARG - Marginalized Kernel Chem. Struct.
#             MINMAX - MinMax Kernel Chem. Struct.
#             SIMCOMP - Graph kernel Chem. Struct.
#             SIDER - Side-effects Similarity Side-effects
#             SPEC - Spectrum Kernel Chem. Struct.
#             TAN - Tanimoto Kernel Chem. Struct.

# Proteins    GIP - Gaussian Interaction Profile Network
#             GO - Gene Ontology Semantic Similarity Func. Annot.
#             MIS-k3m1 - Mismatch kernel
#             (k = 3, m = 1)
#             Sequences
#             MIS-k4m1 - Mismatch kernel
#             (k = 4, m = 1)
#             Sequences
#             MIS-k3m2 - Mismatch kernel
#             (k = 3, m = 2)
#             Sequences
#             MIS-k4m2 - Mismatch kernel
#             (k = 3, m = 2)
#             Sequences
#             PPI - Proximity in protein-protein
#             network
#             Protein-protein
#             Interactions
#             SPEC-k3 - Spectrum kernel (k = 3) Sequences
#             SPEC-k4 - Spectrum kernel (k = 4) Sequences
#             SW - Smith-Waterman aligment score Sequences