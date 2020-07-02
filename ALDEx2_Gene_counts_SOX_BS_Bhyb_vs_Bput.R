  # ALDEx2 on gene abundances of SOX symbionts based on hybrid vs puteoserpentis 2020-01-14
  # Merle Ücker
  
  setwd ("~")
  
  library(ALDEx2)
  library(data.table)
  
  # Import and reformat table
  counts <- fread(file = "Gene_counts_kallisto_SOX_symbionts_Broken_Spur.matrix", header = TRUE, sep = "\t")
  rows <- counts [,1]
  counts <- counts[, c(2:31)]
  counts <- round(counts)
  #counts <- counts [, -c(1:5)] # remove the columns CHROM, POS, START etc. ONLY WITH FEATURECOUNTS DATA!
  
  # Define conditions
  conds <- c(rep("Bhyb", 2), rep("Bput", 2), rep("Bhyb", 1), rep("Bput", 1), rep("Bhyb", 4), rep("Bput", 1), rep("Bhyb", 2), rep("Bput", 3), rep("Bhyb", 2), rep("Bput", 1), rep("Bhyb", 3), rep("Bput", 2), rep("Bhyb", 4), rep("Bput", 2))
  
  # Run ALDEx2
  counts_clr = aldex.clr(counts, conds, mc.samples=128, denom="all", verbose=TRUE, useMC=TRUE)
  counts_kw <- aldex.kw(counts_clr, useMC=TRUE, verbose=TRUE)  
  
  # Identify which values are significant in glm and kw test
  #found.by.all.0.001 <- which(counts_kw$glm.eBH < 0.001 & counts_kw$kw.eBH < 0.001)
  found.by.all.0.05 <- which(counts_kw$glm.eBH < 0.05 & counts_kw$kw.eBH < 0.05)