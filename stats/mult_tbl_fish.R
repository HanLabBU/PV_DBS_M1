# Read in csv file
mod_path <- "/home/pierfier/Dropbox/RKC-HanLab/Pierre PV DBS Project Dropbox/Materials/Plots/Neuronwise/"

# Old file
#df <- read.csv(paste0(mod_path, "modulation_stats.csv"), sep=",", header=TRUE, row.names=1)

# Neuronwise modulation csv without non-modulated neurons included
df <- read.csv(paste0(mod_path, "modulation_stats_keepNonMod0.csv"), sep=",", header=TRUE, row.names=1)

#TODO need to ensure that I am using the correct column names here and the correct rows, because they may have switched

# Printing out the column names
print(colnames(df))
print("")

## Perform fisher's test for transient Vm data
#vm_data <- df[1:(nrow(df)-2), c("Vm.Trans.Act", "Vm.Trans.Sup")]
#
#fish_result <- fisher.test(as.matrix(vm_data))
#print("--- Transient Vm all brain region Fishers ---")
#print(vm_data)
#print(fish_result)
#
##Compare the M1 PV neurons between 140 Hz and 40 Hz
#vm_m1_data <- df[3:4, c("Vm.Trans.Act", "Vm.Trans.Sup")]
#print("--- M1 Vm Transient ---")
#print(vm_m1_data)
#fish_M1_result <- fisher.test(as.matrix(vm_m1_data))
#print(fish_M1_result)
#
##Compare the V1 PV neurons between 140 Hz and 40 Hz
#vm_v1_data <- df[5:6, c("Vm.Trans.Act", "Vm.Trans.Non", "Vm.Trans.Sup")]
#print("--- V1 Vm Transient ---")
#print(vm_v1_data)
#fish_V1_result <- fisher.test(as.matrix(vm_v1_data))
#print(fish_V1_result)
#
## Perform fisher's test for transient firing rate data
#fr_data <- df[3:nrow(df), c("FR.Trans.Act", "FR.Trans.Non", "FR.Trans.Sup")]
#
#fish_result <- fisher.test(as.matrix(fr_data))
#print("--- Transient FR all brain regions Fishers---")
#print(fr_data)
#print(fish_result)
#
## Perform fisher's test for sustained Vm data
#vm_data <- df[3:nrow(df), c("Vm.sus.Act", "Vm.sus.Non", "Vm.sus.Sup")]
#
#fish_result <- fisher.test(as.matrix(vm_data))
#print("---Sustained Vm Fishers---")
#print(vm_data)
#print(fish_result)
#
## Perform fisher's test for sustained firing rate data
#fr_data <- df[3:nrow(df), c("FR.sus.Act", "FR.sus.Non", "FR.sus.Sup")]
#
#fish_result <- fisher.test(as.matrix(fr_data))
#print("---Sustained FR Fishers---")
#print(fr_data)
#print(fish_result)
#

#-- Perform fisher's test for entrained percentages of 140 Hz
rows_140 <- rownames(df)[grep(".*M1.*140|.*V1.*140", rownames(df))]
print(rows_140)
print(rownames(df))
entrain_data <- df[rows_140, c("PLV.Etrain", "PLV.Non.Etrain")]

fish_result <- fisher.test(as.matrix(entrain_data))
print("-- Entrained Fisher's---")
print(entrain_data)
print(fish_result)

#-- Motor cortex entrainment across stimulation frequencies
m1_rows <- rownames(df)[grep("*M1*", rownames(df))]
m1_etrain_data <- df[m1_rows, c("PLV.Etrain", "PLV.Non.Etrain")]

fish_result <- fisher.test(as.matrix(m1_etrain_data))

print("-- M1 Entrained Fisher's --")
print(m1_etrain_data)
print(fish_result)

#-- Perform fisher's test for total Vm activation, Vm suppression, and Vm unchanged across all conditions
# Grab rows for each M1 and V1 condition
cort_rows <- rownames(df)[grep(".*M1.*|.*V1.*", rownames(df))]
print(cort_rows)
#vm_data <- df[cort_rows, c("Vm.Total.Activated", "Vm.Total.Supressed", "Vm.Total.Unchanged")]
vm_data <- df[cort_rows, c("Vm.Total.Activated", "Vm.Total.Supressed")]
#vm_data <- df[cort_rows, c("Vm.Total.Supressed")] # This does not work
fish_result <- fisher.test(as.matrix(vm_data), workspace = 2e8)
print("--Vm Total All Conditions Fisher's--")
print(vm_data)
print(fish_result)
