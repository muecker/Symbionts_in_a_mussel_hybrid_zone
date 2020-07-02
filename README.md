# Symbionts in a mussel hybrid zone
The purpose of this repository is to illustrate the analyses performed in our study about sulfur-oxidising symbionts of hybrid and parental Bathymodiolus mussels in Broken Spur. The scripts should help interested readers to better understand the methods used. Please be aware that they are not tools intended for further usage.

## Analysis of average nucleotide identity of Bathymodiolus SOX symbionts from Broken Spur
_Visualise ANI values of SOX symbionts at Broken Spur and test for correlation with host genetics and sampling year_ \
Script: `Heatmap+Mantel_ANI_SOX_BS.R` \
Data needed: Average_nucleotide_identity_SOX_symbionts_Broken_Spur.csv, Pairwise_distances_host_SNPs_ANI.txt

## Analysis of FST values of Bathymodiolus SOX symbionts from Broken Spur
_Visualise FST values of SOX symbionts at Broken Spur and test for correlation with host genetics and sampling year_ \
Script: `Heatmap+Mantel_FST_SOX_BS.R` \
Data needed: Pairwise_mean_FST_SOX_symbionts_Broken_Spur.csv, Pairwise_distances_host_SNPs_FST.txt

## Analysis of Bathymodiolus SOX symbiont gene counts from Broken Spur
_Test for difference in gene abundances between symbionts of hybrids and parental mussels with ALDEx2_ \
Script: `ALDEx2_Gene_counts_SOX_BS_Bhyb_vs_Bput.R` \
Data needed: Gene_counts_kallisto_SOX_symbionts_Broken_Spur.matrix

## Analysis of Bathymodiolus SOX symbiont FST vs. geography at the NMAR (addition to phylogenomic analysis)
_Plot and test for correlation between symbiont FST and geographic distances_ \
Script: `Mantel+Plot_FST_vs_geographic_distances_SOX_symbionts_NMAR.Rmd` \
Data needed: Phylogenomics_NMAR.fasta, Groups_NMAR.csv, Pop_NMAR.csv, Geographic_distances_NMAR.csv

## Redundancy analysis of Bathymodiolus SOX symbiont allele frequencies at the NMAR
_Analyse the variation in symbiont allele frequencies that can be explained by either spatial information, host species, rock type or depth_ \
Script: `Redundancy_analysis_allele_frequencies.Rmd` \
Data needed: Allele_frequencies_SOX_symbionts_NMAR.csv, Metadata_for_RDA.csv

## Analysis of Bathymodiolus SOX symbiont per gene FST values
_Find genes that are more differentiated between symbionts of hybrids and parental mussels than within either of the groups_ \
Script: `Per_gene_FST_analysis.Rmd` \
Data needed: Per_gene_FST_SOX_symbionts_Broken_Spur.zip
