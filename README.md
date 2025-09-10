# MH_orchid_mycorrhizae_phyloseq
Code for processing ITS2 metabarcoding data from Corallorhiza spp. and Cephalanthera austiniae
The following code processes ITS2 data on orchid mycorrhizas.

- Read Processing, clustering, denoising, etc. was conducted with PIPITS using SINTAX and the UNITE fungal ITS2 reference database:
- Gweon, H. S., Oliver, A., Taylor, J., Booth, T., Gibbs, M., Read, D. S., Griffiths, R. I., & Schonrogge, K. (2015). PIPITS: an automated pipeline for analyses of fungal internal transcribed spacer sequences from the Illumina sequencing platform. Methods in Ecology and Evolution, 6(8), 973–980. https://doi.org/10.1111/2041-210X.12399
- Imports OTU table (biom) and metadata with phyloseq
- Uses TMM normalization with edgeR instead of rarefaction
- Filter out all but Kingdom = Fungi
- Runs decontam package to remove contaminants
- Removes extremely rare OTUs that might be errors
- Standardizes taxonomy/rank
- Filters all potential orchid, ericoid, ecto-, and arbuscular mycorrhizal taxa by family
- Allows focal species (of orchid) selection, e.g. Corallorhiza striata, Cephalanthera austiniae, etc.
- Assigns Trophic Mode to each OTU with FUNGuildR
- Plots relative abundance with ggplot2 and plotLy interactively
