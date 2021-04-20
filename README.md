# scRNAmod
Registration of nanowell scanning images and SeqFISH as described in the publication of "Systematic detection of m6A-modified transcripts at single-molecule and single-cell resolution" 

# Dependencies
This software requires
1. R
2. Python
3. ImageJ/Fiji

# Scripts
1. NanoWell_FOV_Alignment.R - convert the X/Y coordinates of each nanowell along with the number of loaded cells from the fluorescence image into a cell occupancy matrix (COM)
2. Register_NanoWell_FOV_Info.R - register an RNA density matrix (RDM) of transcript abundance and m6A modification levels derived from TIRF images with COM
3. SeqFISH_Gene_Table.R - load SeqFISH data into a gene expression matrix (GEM)
4. Merge_Quant.R - add GEM to the registered COM subset
5. Merge_CellType.R - add cell phenotypes to the registered COM subset
6. tSNE.R - unsupervised mapping of cell types from SeqFISH gene expressions
7. SeqFISH_alignment.py - align images and count the number of SeqFISH signals
8. Nanowell_Intensity.ijm - single cell phenotyping in multicolor nanowell scanning images

# Updates
updated on 04/02/2021

# Contact
Kyung Lock Kim (klkim@mgh.harvard.edu)
