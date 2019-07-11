# Calculate mutual information for protein sequence alignment of Dpr-DIP complexes

Our aim is to calculate mutual information for interface contacts of Dpr-DIP protein complex. 
In particular, we want to compare top mutual information for interface contacts between cognate 
and non-cognate complexes.

The calculation of mutual information involves the following steps: calculate mutual information 
between each residue of Dpr and each residue of DIP. From this set, we filter out the residue 
indices that are present at the interface for each protein and sort the values of mutual 
information for each interface contacts. The next step is to map the indices of residues in 
each pdb file to indices of residues in the sequence alignment file. In the last step we save the 
pdb and MSA indices and their corresponding mutual information for each protein complex. 

For visualization purposes, we form an union set of top 10 mutual information for cognate and 
non-cognate complexes and visualize in VMD. 

We are currently looking at ways to distinguish the pairwise interaction details of cognate and non-cognate set.

TODO: 1) Acquire larger data sets for DIP/DPR interactions
      2) Bootstrap the original data sets to resample for better statistics
      3) Implement expectation value based statistics to compare cognate / non-cognate interactions