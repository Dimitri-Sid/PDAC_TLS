# - - - - - - - - - - - - 
# Process CODA output
# Author: Dimitri Sidiropoulos
# - - - - - - - - - - - - 

source("./0_Visium_config_functions.R")

CODA_output<- "./rawData/CODA_output"
samples <- c("J1568_C1","J1568_C2","J1568_C3","J1568_C4_manual_alignment","J1568_C5","J1568_C6","J1568_C7_manual_alignment","J1568_C8","J1568_C9","J1568_C10","J1568_C11","J1568_C12")
samplespaths <- paste0(CODA_output,"/",samples)

## Combine all CODA outputs into single objects
coda_cell_props <- data.frame()
for (i in 1:12){
  coda <- data.frame(read_excel(paste0(CODA_output, "/",samples[i],"/tissue_positions_cellular_compositions.xlsx")))
  rownames(coda) <- paste0(samples[i],"_",coda$...1)
  coda$sample <- samples[i]
  coda_cell_props <- rbind(coda_cell_props,coda)
}

coda_cell_count <- data.frame()
for (i in 1:12){
  coda <- data.frame(read_excel(paste0(CODA_output, "/",samples[i],"/tissue_positions_cellular_count.xlsx")))
  rownames(coda) <- paste0(samples[i],"_",coda$...1)
  coda$sample <- samples[i]
  coda_cell_count <- rbind(coda_cell_count,coda)
}

coda_tissue_props <- data.frame()
for (i in 1:12){
  coda <- data.frame(read_excel(paste0(CODA_output, "/",samples[i],"/tissue_positions_tissue_compositions.xlsx")))
  rownames(coda) <- paste0(samples[i],"_",coda$...1)
  coda$sample <- samples[i]
  coda_tissue_props <- rbind(coda_tissue_props,coda)
}

write.table(coda_cell_props, file = "./dataObjects/coda_cell_props_all_samples.txt")
write.table(coda_cell_count, file = "coda_cell_count_all_samples.txt")
write.table(coda_tissue_props, file = "coda_tissue_props_all_samples.txt")
