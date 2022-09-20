
The Woolly marmot dataset (simulated data) -- Thibaut Capblancq

In woolly marmots, three quantitative traits are particularily involved in local adaptation - wool length, wool color and production of lanolin (wool wax). Each one of these traits is coded by a relatively small number of genes and the precise environmental variables driving selection are still unknown but should be among the 10 following variables: 

- BWS: Bad Weather all Summer
- TSP: Tons of Skiing Potential
- SCS: Super Cold Summer
- WR: Whiteness of the Rocks
- GG: Greenness of the Grass
- NTS: Number of Tourists in Summer
- LSS: Lots of Snow during Summer 
- PE: Presence of Eagles
- SCW: Super Cold Winter
- DR: Density of Rocks

Ten individuals were sampled at 61 source localities. The woolly_marmot_data.txt file contains both meta data and genotypes for 610 individual samples, with individual samples in rows and environmental variables and genetic loci in columns. Individuals were genotyped at 1000 diploid loci, which are encoded in genepop format: 0 for an ancestral allele homozygote, 1 for an heterozygote and 2 for a derived allele homozygote. No missing data.

Note: the *woollymarmot.baypass.geno* file is the allele count data file in *BayPass* format (i.e., it contains allele count for the reference and alternate alleles in each of the 61 populations for the 1,000 SNPs). The order of the SNP and population is the same as in the original file. Likewise, the *woollymarmot.baypass.cov* is the population covariable format in *BayPass* format (i.e., it contains the 10 environmental covariable values for each of the 61 populations).

The Woolly_marmot_reintroduction.txt file contains environmental information for the three sites selected for reintroducing the Woolly marmot.
