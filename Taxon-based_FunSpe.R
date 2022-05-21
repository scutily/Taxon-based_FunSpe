# A ecosystem function that utilizing Glucose-1-phosphate 
EF <- read.csv("EF_G1P.csv")  

# Community omposition (ASV table) of prokaryotes.
composition_ASV <- read.csv("composition_ASV.csv")
row.names(composition_ASV) <- composition_ASV[,1]
composition_ASV <- as.matrix(composition_ASV[,-1])


###############################################
# Perform the Msa-enet based on modified code #
###############################################
# Please run the code "msaenet.new.R" before running the code below

# install.packages("doParallel")
library("doParallel")
registerDoParallel(detectCores()) ## Speed up the computation

a <- c(1:9)/10   # Set alpha values from 0.1-0.9 for sparse regression (i.e., pure elastic net)


msaenet_list <- msaenet.new(composition_ASV,EF[,2],"gaussian","ridge" # Initial step = ridge
                           ,alpha=a,  nfolds=nrow(composition_ASV), # i.e., Leave-one-out cross-validation
                           rule="lambda.min",seed=1010,
                           nsteps = 10L, # 10 steps for variables selection
                           tune.nsteps = "ebic",parallel = TRUE) 
msaenet_mse <- msaenet_list$min.errors # Extract the mse for each step
msaenet_minmse_coefficients<- msaenet_list$beta.list[[which.min(msaenet_mse)]] # Select the step depending on the minimum mse
  
abs_faci_ASV <- length(which(msaenet_minmse_coefficients>0))  # Absolute number of facilitative ASV
abs_repr_ASV <- length(which(msaenet_minmse_coefficients<0))  # Absolute number of repressive ASV
prop_faci_ASV <- abs_faci_ASV/(abs_faci_ASV+abs_repr_ASV) # Proportion of facilitative ASV

specificity_taxon <- data.frame(abs_faci_ASV, abs_repr_ASV, prop_faci_ASV)
specificity_taxon
