library(readr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(caret)

get_control <- function(data_path="JuliaSero 020318 1h drug treatments/single cells by condition (3 replicates per file)",
                        corr_threshold=0.99){
  idx_removed = c(1:10, 14, 24, 30, 34, 38, 42, 
                  55:57,
                  61, 70:72)
  
  Control_Nuclei <- read.csv(sprintf("%s/Control_Nuclei.csv", data_path))
  
  df_Control <- Control_Nuclei %>% .[-idx_removed] %>%
    na.omit() %>% as.data.frame() %>%
    .[-c(7, 12, 20, 21, 29, 33)]
    # select_col_by_cor(df=., corr_threshold)
  
  return(df_Control)
}

get_treatment <- function(data_path="JuliaSero 020318 1h drug treatments/single cells by condition (3 replicates per file)",
                          ctl_vars, treatment_method="VS4718_Nuclei"){
  Treatment_Nuclei <- read.csv(sprintf("%s/%s.csv", data_path,
                                       treatment_method))
  df_Treatment <- Treatment_Nuclei %>%
    dplyr::select(all_of(ctl_vars)) %>%
    na.omit() %>% as.data.frame()
  
  return(df_Treatment)
}


# support functions --------------------------------
library(dendextend)
library(heatmaply)
library(tibble)



select_col_by_cor <- function(df, threshold=0.9){
  cor_abs_ctl = df %>%
    # select(-c(treatment)) %>%
    cor() %>% abs()
  map_c = cor_abs_ctl%>%
    heatmap(keep.dendro=TRUE)
  dend1 = map_c$Rowv
  # source("./week5 highly correlated rm.R")
  
  corr_threshold = threshold
  k = get_min_nclust(df, corr_threshold, dend_ = dend1)
  clustered_group_ctl = cut_dend_get_group(k, dend = dend1)
  # print(clustered_group_ctl)
  df %>% dplyr::select(clustered_group_ctl$repre) %>% return()
}


# get clustered groups
cut_dend_get_group <- function(cluster_num=34, 
                               dend=dend1){
  
  result1 = cutree(dend, k=cluster_num)
  result1 <- result1 %>% as.data.frame() 
  colnames(result1) <- c("hc_group")
  result1 = result1%>%
    tibble::rownames_to_column(., "Variables")
  # print(result1)
  cluster_table1 = result1 %>% group_by(hc_group) %>% 
    mutate(vars_by_hc = paste0(Variables, collapse = "; "),
           repre = Variables[1]) %>% 
    dplyr::select(hc_group, vars_by_hc, repre) %>%
    unique()
  # print(cluster_table1)
  
  return(cluster_table1)
  
}



get_min_nclust = function(df, p, dend_=dend1, is_print_=FALSE){
  for(k in 25:50){
    if(check_corr(df, k, threshold = p, dend=dend_, is_print=is_print_)){
      if(is_print_){
        print(sprintf("Successful at cluster number = %i", k))
      }
      return(k)
    }
  }
  return(-1)
}

check_corr <- function(df, cluster_num, threshold=0.9, 
                       dend=dend1,
                       is_print=FALSE){
  
  result1 = cutree(dend, k=cluster_num)
  result1 <- result1 %>% as.data.frame() 
  colnames(result1) <- c("hc_group")
  result1 = result1%>%
    tibble::rownames_to_column(., "Variables")
  
  cluster_table1 = result1 %>% group_by(hc_group) %>% 
    mutate(vars_by_hc = paste0(Variables, collapse = ";")) %>% 
    dplyr::select(hc_group, vars_by_hc) %>%
    unique() %>%
    .$vars_by_hc
  
  for(i in 1:cluster_num){
    variables = str_split(cluster_table1[i], pattern=";")[[1]]
    if(length(variables)>1){
      corr_i = df %>%   
        dplyr::select(variables) %>%
        cor(method = "pearson") %>% abs()
      if(all(corr_i>=threshold)==FALSE){
        if(is_print){
          print(sprintf("cluster number %s not successful for group %i: %s", 
                        cluster_num, i, cluster_table1[i]))
        }
        
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


get_category <- function(var_name){
  if(startsWith(var_name, "AreaShape")){
    return("NucMorph")
  } else if (startsWith(var_name, "Intensity")){
    return("NucProLife")
  } else if(var_contains_char(c("Math_Cell" ,"Area"), var_name)){
    return("CellMorph")
  } else if (var_contains_char(c("Density" ,"RadialMean", "Intensity"), var_name)){
    return("CellProLife")
  } else{
    return("Others")
  }
}

# colnames(x) %>% get_category()

var_contains_char <- function(char_list, variable_name){
  char_list %>% 
    sapply(., function(c){grepl(c, variable_name, fixed = TRUE)}) %>% 
    any() %>% 
    return()
}

# var_contains_char(c("a" ,"b", "c"), "ab")
