# Get pairwise genetic distance of mussel host from Fluidigm Genotyping data 2020-02-05
# Merle Ücker

# Load libraries
library(ape)
library(adegenet)

# Load genotyping data and get genetic distances
snp_str <- read.structure(file = "Host_genotyping.str")
snp_str_df <- genind2df(snp_str)
distances_str_ape <- as.matrix(dist.gene(snp_str_df, method="percentage"))
distances_str_ape_BS <- as.matrix(distances_str_ape[111:153, 111:153])

# Load names
names <- read.table("names_for_host_distances")
rownames(distances_str_ape_BS) <- names[,1]
colnames(distances_str_ape_BS) <- names[,1]

# Reorder genetic distances and save them for further analysis with symbiont ANI values
col.order <- c("2424_D_Bhyb","2424_E_Bhyb","2424_H_Bput","2424_N_Bhyb","3386_AB_Bput","3386_AD_Bhyb","3386_AE_Bput","3386_A_Bput","3386_AG_Bhyb","3386_AH_Bhyb","3386_AJ_Bhyb","3386_AK_Bput","3386_AL_Bhyb","3386_C_Bhyb","3386_D_Bput","3386_E_Bput","3386_F_Bput","3386_G_Bhyb","3386_H_Bhyb","3386_I_Bput","3386_J_Bhyb","3386_M_Bhyb","3386_O_Bput","3386_S_Bput","3386_T_Bhyb","3386_U_Bhyb","3386_W_Bhyb","3386_X_Bhyb","3386_Y_Bput","3386_Z_Bput")
row.order <- c("2424_D_Bhyb","2424_E_Bhyb","2424_H_Bput","2424_N_Bhyb","3386_AB_Bput","3386_AD_Bhyb","3386_AE_Bput","3386_A_Bput","3386_AG_Bhyb","3386_AH_Bhyb","3386_AJ_Bhyb","3386_AK_Bput","3386_AL_Bhyb","3386_C_Bhyb","3386_D_Bput","3386_E_Bput","3386_F_Bput","3386_G_Bhyb","3386_H_Bhyb","3386_I_Bput","3386_J_Bhyb","3386_M_Bhyb","3386_O_Bput","3386_S_Bput","3386_T_Bhyb","3386_U_Bhyb","3386_W_Bhyb","3386_X_Bhyb","3386_Y_Bput","3386_Z_Bput")
distances_str_ape_BS_reordered <- as.matrix(distances_str_ape_BS[row.order,col.order])
write.table(distances_str_ape_BS_reordered, file = "pairwise_distances_host_SNPs_ANI.txt")

# Reorder genetic distances and save them for further analysis with symbiont FST values
fst_col.order <- c("2424_E_Bhyb","2424_N_Bhyb","3386_A_Bput","3386_AB_Bput","3386_AD_Bhyb","3386_AE_Bput","3386_AG_Bhyb","3386_AH_Bhyb","3386_AI_Bhyb","3386_AJ_Bhyb","3386_AK_Bput","3386_AL_Bhyb","3386_C_Bhyb","3386_D_Bput","3386_E_Bput","3386_F_Bput","3386_G_Bhyb","3386_H_Bhyb","3386_I_Bput","3386_J_Bhyb","3386_M_Bhyb","3386_N_Bhyb","3386_O_Bput","3386_S_Bput","3386_T_Bhyb","3386_U_Bhyb","3386_W_Bhyb","3386_X_Bhyb","3386_Y_Bput","3386_Z_Bput")
fst_row.order <- c("2424_E_Bhyb","2424_N_Bhyb","3386_A_Bput","3386_AB_Bput","3386_AD_Bhyb","3386_AE_Bput","3386_AG_Bhyb","3386_AH_Bhyb","3386_AI_Bhyb","3386_AJ_Bhyb","3386_AK_Bput","3386_AL_Bhyb","3386_C_Bhyb","3386_D_Bput","3386_E_Bput","3386_F_Bput","3386_G_Bhyb","3386_H_Bhyb","3386_I_Bput","3386_J_Bhyb","3386_M_Bhyb","3386_N_Bhyb","3386_O_Bput","3386_S_Bput","3386_T_Bhyb","3386_U_Bhyb","3386_W_Bhyb","3386_X_Bhyb","3386_Y_Bput","3386_Z_Bput")
distances_str_ape_BS_reordered_fst <- as.matrix(distances_str_ape_BS[fst_row.order,fst_col.order])
write.table(distances_str_ape_BS_reordered_fst, file = "pairwise_distances_host_SNPs_FST.txt")