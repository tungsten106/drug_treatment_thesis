#week8 hparam tuning
library(readr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(caret)
library(randomForest)
library(tictoc)
library(doParallel)
library(foreach)

registerDoParallel(cores=detectCores())


source("./week8 functions.R")
source("./week8 preprocessing.R")

data_path = "JuliaSero 020318 1h drug treatments/single cells by condition (3 replicates per file)"
treatments_list = list.files(path=data_path, pattern=NULL, all.files=FALSE,
                             full.names=FALSE) %>% tools::file_path_sans_ext() %>%
  .[. != "Control_Nuclei"]

df_Control = get_control(data_path=data_path)

# tuning start ------------------------------------------------------------

new_folder1 = sprintf("hparams_wk8_var44")

dir.create(sprintf("./%s", new_folder1))

for(treatment_method in c("Control_Nuclei",treatments_list)){
  tic(treatment_method)
  df_Treatment = get_treatment(data_path=data_path, 
                               ctl_vars=colnames(df_Control),
                               treatment_method=treatment_method)
  
  # get predictor variables and response variables
  x1 <- subset(df_Treatment[, ], select = c(-Math_logNucByRingYAP))
  y1 <- (subset(df_Treatment[, ], select = Math_logNucByRingYAP))
  y1 <- y1[,1]
  
  # tune hparam by oob-mse
  rf_select_result_1_parall = 
    get_param_parallel(x1, y1, 
                       SUBSET_SIZE = 5000, # tune with 5000 subset
                       ntrees = c(100, 200,300,400,500),
                       mtrys=seq(20, ncol(x1), by=1),
                       nodesizes=seq(1,5))
  
  
  param_tr = rf_select_result_1_parall$mse_grid %>% arrange(oob_mses) %>% .[1,]
  
  # save mse grid just in case for plotting
  write.csv(rf_select_result_1_parall$mse_grid, 
            sprintf("./%s/mse_%s.csv", new_folder1, treatment_method), 
            row.names = FALSE) 
  
  # write.csv(param_tr, sprintf("./%s/hp_%s.csv", new_folder1, treatment_method),
  #           row.names = FALSE)
  
  toc() #  with parallel
}

# conbine them together ---------------------

tuned_hparam_list = data.frame()
# tuned_hparam_list %>% rbind(., best_hparam)

for(treatment_method in c("Control_Nuclei", treatments_list)){
  mse_grid = read.csv(sprintf("./hparams/mse_%s.csv", treatment_method))
  
  best_hparam = mse_grid %>% arrange(oob_mses) %>% .[1,]
  best_hparam$treatment = treatment_method
  
  tuned_hparam_list <- tuned_hparam_list %>% rbind(., best_hparam)
}

write.csv(tuned_hparam_list, sprintf("./%s/hparam_all.csv", new_folder2), 
          row.names = FALSE)

# train to get CATE -------------------------------
new_folder2 = "CATE_wk8_var44"

hparams = read.csv(sprintf("./%s/hparam_all.csv", new_folder2))
phi_0 = hparams %>% filter(treatment == "Control_Nuclei")

dir.create(sprintf("./%s", new_folder2))

for(treatment_method in treatments_list){
  df_Treatment = get_treatment(data_path=data_path, 
                               ctl_vars=colnames(df_Control),
                               treatment_method=treatment_method)
  x1 <- subset(df_Treatment[, ], select = c(-Math_logNucByRingYAP))
  y1 <- (subset(df_Treatment[, ], select = Math_logNucByRingYAP))
  y1 <- y1[,1]
  
  
  # param_1 = read.csv(sprintf("./hparams/hp_%s.csv", treatment_method))
  phi_1 = hparams %>% filter(treatment == treatment_method)
  
  tic(treatment_method)
  CATE_tr = T_learner(x0, y0, x1, y1,
                      hparam0 = phi_0,
                      hparam1 = phi_1)
  toc() 
  
  write.csv(CATE_tr,sprintf("./%s/CATE_%s.csv", new_folder2, treatment_method), 
            row.names = FALSE)
}

# plot CATE --------------

CATE_all = data.frame()
for(treatment_method in treatments_list){
  cate_tmp = read.csv(sprintf("./%s/CATE_%s.csv", new_folder2, treatment_method))
  meltData = melt(cate_tmp)
  meltData$treatment = treatment_method
  CATE_all = rbind(CATE_all, meltData)
}

pdf(sprintf("./%s/CATE_all_plot.csv", new_folder2))
CATE_all %>% 
  ggplot(aes(factor(variable), value))  +
  geom_boxplot() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=0.8) +
  facet_wrap(~treatment)
dev.off()