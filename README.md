# the-cut
Analyses an amino acid sequence and returns where a specific protease is likely to cut. The program performs a frequency analysis of the amino acid residues preceeding and following the cut position (-4 to +4) based on all known substrates (using data from the Merops database). The amino acid sequence is analyzed in 8 amino acid sequence windows and a numerical probability value is assigned to each window and the windows are returned in the order of highest probability.

We use the meropsrefs.sql table from https://www.ebi.ac.uk/merops/
