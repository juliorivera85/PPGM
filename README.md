# PPGM
Basics on reconstructing ancestral suitable habitats

For detail read:

Rivera, Lawing, and Martins (2020) Reconstructing historical shifts in suitable habitat of Sceloporus lineages using phylogenetic niche modeling. Journal of Biogeography.

This is an R analysis that allows you to reconstruct suitable habitat for a species in North and Central America from 1 to 20 million years ago at 1 MY increments.

The main code is called script_PPGM_JR and there are various smaller piecies of code that need to be sourced in a file called source_function

I've also included code to run an ancestral reconstruction analysis of the climate variables along the phylogeny within the PPGM code. Then, you can find the phylogenetic signal of that climate variable you just reconstrcuted.

Last, there is code to  run an LTT analysis that will help you understand how speciation rates vary through time. 

I've inluded the paleocliamte values in a Rdata file called paleoclimate. You don't have to modify this. 
