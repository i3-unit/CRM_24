CORES <- parallel::detectCores()/2

library(entropy)
library(dplyr, warn.conflicts = FALSE)
library(data.table)
library(reshape2)
library(stringdist)
library(stringr)
library(ggraph)
library(tidyr)
library(RColorBrewer)
library(igraph)
library(stringdist)
library(gdata)
library(grid)
library(ggsci)
library(iNEXT)
library(ggplot2)
library(Rarity)
library(ggthemes)
library(doParallel)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(ggpattern)
library(ComplexHeatmap)
library(readxl)
select = dplyr::select


read_files <- function(dt.sample, nt=FALSE) {
  if(nt==TRUE){
    x <- fread(dt.sample, select = c("cloneId", "cloneCount","nSeqCDR3", "aaSeqCDR3", "allVHitsWithScore", "allJHitsWithScore"))
  } else {
  x <- fread(dt.sample, select = c("cloneId", "cloneCount", "aaSeqCDR3", "allVHitsWithScore", "allJHitsWithScore"))
  }
}


named_group_split <- function(.tbl, ..., sep = sep) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = sep)))
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

data_modelling <- function(data, row, col, value) {
  data = setDT(data) %>% dcast.data.table(get(row) ~ get(col), value.var = value, fill=0, fun = sum)
  data = as.data.frame(data)
  rownames(data) = data$row
  data = data %>% dplyr::select(-row)
  data = as.matrix(data)
}

StatCentSeg <- ggplot2::ggproto("StatCentSeg", Stat,
                                compute_group = function(data, scales, params,
                                                         cfun=median) {
                                  data$xend <- cfun(data$x)
                                  data$yend <- cfun(data$y)
                                  return(data)
                                },
                                required_aes = c("x", "y")
)
stat_centseg <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA, 
                         inherit.aes = TRUE, cfun=median, ...) {
  layer(
    stat = StatCentSeg, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, cfun = cfun, ...)
  )
}


##### First data processing 
prefilter_df <- function (data, var, nt) {
 
   message_parallel <- function(...){
    system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
  }
  message_parallel(unique(data$filename))
  
   data <- data %>%
    dplyr::select(filename, nSeqCDR3, cloneCount, allVHitsWithScore, allJHitsWithScore, aaSeqCDR3) %>%
    setnames(old = c("aaSeqCDR3"), new = c("cdr3aa"))
  
  data$cloneCount <- as.numeric(data$cloneCount)
  #data <- data %>% filter(cloneCount >= var$ana$cut_off_cloneCount)
  
  # Remove the non-productive and the longest sequences
  data <- data %>% 
    filter(!str_detect(cdr3aa, "[*]")) %>%
    filter(!str_detect(cdr3aa, "[_]")) %>%
    mutate(length_cdr3aa = nchar(cdr3aa)) %>% 
    filter(length_cdr3aa > mean(length_cdr3aa)-8 & length_cdr3aa < mean(length_cdr3aa)+8) %>%
    filter(!cdr3aa == "YAEQFF")
  
  
  # Add column TRV and TRJ, remove .txt in filename
  data <- data %>% 
    separate(allVHitsWithScore, into = "TRV", sep = "[*]") %>%
    separate(allJHitsWithScore, into = "TRJ", sep = "[*]") #%>%
  #separate(filename, into = "filename", sep = "[.]")
  
  data <- as.data.table(data)
  
  # Create column Chain to identify TRA or TRB
  data <- data %>% 
    filter(!str_detect(TRV, 'TRD')) %>%   # we look in the column TRV only the raw who don't contain 'TRD' 
    mutate(chain = if_else(str_detect(TRV, 'TRA') == "TRUE", "TRA",
                           if_else(str_detect(TRV, 'TRB') == "TRUE", "TRB", "others")))
  
  # Remove chain we don't want, like : IGKV4-1, IGKJ1, TRGV5, etc.
  data <- data[data$chain != "others",]
  
  # Create Clonotype column
  if(nt==TRUE){
    data$clonotypes <- paste(data$TRV, data$nSeqCDR3,data$cdr3aa, data$TRJ)
  } else {
    data$clonotypes <- paste(data$TRV, data$cdr3aa, data$TRJ)
  }
  
  data$VJ <- paste(data$TRV, data$TRJ)
  
  
  # Calculate the frequency by sample and chain
  data <- data %>%
    group_by(clonotypes, filename, chain) %>%
    mutate(cloneCount = sum(cloneCount)) %>%
    unique()
  
  data <- data %>%
    group_by(filename, chain) %>%
    mutate(freq_clo = cloneCount / sum(cloneCount) * 100) %>%
    ungroup()
  
  # Calculate the number of clonotypes by filename
  data <- data %>%
    dplyr::group_by(filename, chain) %>%
    mutate(nb_clonotypes = length(unique(clonotypes))) %>%
    mutate(nb_sequences = sum(cloneCount)) %>%
    mutate(sequences_per_clonotypes = round(nb_sequences/nb_clonotypes, 3))
  

  print(paste0("data row number BEFORE merge with metadata : ", nrow(data)))
  data <- merge(data, metadata_2 <- unique(metadata[,-c("chain")]), by.x = "filename", by.y = "filename", all.x = T)
  print(paste0("data row number AFTER merge with metadata : ", nrow(data)))
  #data$sample_id <- paste(data$filename, data$chain, sep = "/")
  #tmp_files(paste0(message, "A prefilter has been made on the data"), var)

  return(data)
}

theme_Publication <- function(base_size=12, base_family="sans") {

  (theme_light(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "black"),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            #axis.line.x = element_line(colour="black"),
            #axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="black",fill="#f0f0f0"),
            strip.text = element_text(face="bold", colour = "black")
    ))
}


dkast <- function(data, gene, value) {
  #data = data %>% filter(chain == get(chain)) %>% filter(cell_subset_acronym == get(cell_subset))
  data = setDT(data) %>% dcast.data.table(get(gene) ~ group2, value.var = value, fill=0, fun = sum)
  data = data %>% dplyr::select(-gene)
  data = as.matrix(data)
}



Compute_MH <- function(MPC, method, out="matrix") {
  print(colSums(MPC))
  #MPC <- MPC %>% select_if(colSums(.) > 10)
  if (method == "ji") {
    TRA_MH <- as.matrix(divo::ji(MPC))
  }
  if (method == "mh") {
    TRA_MH <- as.matrix(divo::mh(MPC))
  }
  
  TRA_MH <- do.call(rbind.data.frame, TRA_MH)
  TRA_MH <- TRA_MH[1:ncol(TRA_MH),]
  gdata::upperTriangle(TRA_MH, diag = T, byrow=T) <- gdata::lowerTriangle(TRA_MH, diag = T)
  
  if(out== "matrix"){
    return(TRA_MH)
  } else {
  TRA_MH_df = melt(TRA_MH)
  TRA_MH_df = TRA_MH_df %>% filter(value < 1)
  return(TRA_MH_df)
  }
}









###########################################
###########################################
###########################################
###########################################
###########################################
###########################################

CDR3aa_group_function <- function(X_CL) {
  if (var_cell_subset == "nTregs") {
    
    A <- X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 > 0, indiv_3 > 0, indiv_4 > 0, indiv_5 > 0, indiv_6 > 0, indiv_7 > 0, indiv_8 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "all_indiv")
    colnames(A) <- c("cdr3aa", "group")
    
    A2 <- X_CL %>% filter(pool_indiv == 0, indiv_1 >= 0, indiv_2 >= 0, indiv_3 >= 0, indiv_4 >= 0, indiv_5 >= 0, indiv_6 >= 0, indiv_7 >= 0, indiv_8 >= 0, indiv_9 >= 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "At_least_2_indiv")
    colnames(A2) <- c("cdr3aa", "group")
    
    B <-  X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(B) <- c("cdr3aa", "group")
    
    C <- X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 > 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(C) <- c("cdr3aa", "group")
    
    D <- X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 > 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(D) <- c("cdr3aa", "group")
    
    E <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 > 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(E) <- c("cdr3aa", "group")
    
    G <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 > 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(G) <- c("cdr3aa", "group")
    
    H <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 > 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H) <- c("cdr3aa", "group")
    
    H2 <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 > 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H2) <- c("cdr3aa", "group")
    
    H3 <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 > 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H3) <- c("cdr3aa", "group")
    
    H4 <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H4) <- c("cdr3aa", "group")
    
    I <- X_CL %>% filter(pool_indiv > 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "pool_private")
    colnames(I) <- c("cdr3aa", "group")
    
    J <-  X_CL %>% filter(pool_indiv > 0, indiv_1 > 0, indiv_2 > 0, indiv_3> 0, indiv_4 > 0, indiv_5 > 0, indiv_6 > 0, indiv_7 > 0, indiv_8 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "Shared_by_all")
    colnames(J) <- c("cdr3aa", "group")
    
    Indiv <- rbind(B, C, D, E, G, H, H2, H3, H4)
    A2 = A2 %>% filter(!cdr3aa %in% Indiv$cdr3aa)
    A2 = A2 %>% filter(!cdr3aa %in% A$cdr3aa)
    
    CDR3aa_group = rbind(A, A2, Indiv, I, J)
    colnames(CDR3aa_group) <- c("cdr3aa", "group")
    return(CDR3aa_group)
  }
  
  
  
  if (var_cell_subset == "amTregs") {
    
    A <- X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 > 0, indiv_3 > 0, indiv_4 > 0, indiv_5 > 0, indiv_7 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "all_indiv")
    colnames(A) <- c("cdr3aa", "group")
    
    A2 <- X_CL %>% filter(pool_indiv == 0, indiv_1 >= 0, indiv_2 >= 0, indiv_3 >= 0, indiv_4 >= 0, indiv_5 >= 0, indiv_7 >= 0, indiv_9 >= 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "At_least_2_indiv")
    colnames(A2) <- c("cdr3aa", "group")
    
    B <-  X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(B) <- c("cdr3aa", "group")
    
    C <- X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 > 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(C) <- c("cdr3aa", "group")
    
    D <- X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 > 0, indiv_4 == 0, indiv_5 == 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(D) <- c("cdr3aa", "group")
    
    E <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 > 0, indiv_5 == 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(E) <- c("cdr3aa", "group")
    
    G <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 > 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(G) <- c("cdr3aa", "group")
    
    H2 <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_7 > 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H2) <- c("cdr3aa", "group")
    
    H4 <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_7 == 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H4) <- c("cdr3aa", "group")
    
    I <- X_CL %>% filter(pool_indiv > 0, indiv_1 == 0, indiv_2 == 0, indiv_3 == 0, indiv_4 == 0, indiv_5 == 0, indiv_7 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "pool_private")
    colnames(I) <- c("cdr3aa", "group")
    
    J <-  X_CL %>% filter(pool_indiv > 0, indiv_1 > 0, indiv_2 > 0, indiv_3> 0, indiv_4 > 0, indiv_5 > 0, indiv_7 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "Shared_by_all")
    colnames(J) <- c("cdr3aa", "group")
    
    Indiv <- rbind(B, C, D, E, G, H2, H4)
    A2 = A2 %>% filter(!cdr3aa %in% Indiv$cdr3aa)
    A2 = A2 %>% filter(!cdr3aa %in% A$cdr3aa)
    
    CDR3aa_group = rbind(A, A2, Indiv, I, J)
    colnames(CDR3aa_group) <- c("cdr3aa", "group")
    return(CDR3aa_group)
  }
  
  
  
  
  if (var_cell_subset == "CD8") {
    
    A <- X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 > 0, indiv_6 > 0, indiv_7 > 0, indiv_8 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "all_indiv")
    colnames(A) <- c("cdr3aa", "group")
    
    A2 <- X_CL %>% filter(pool_indiv == 0, indiv_1 >= 0, indiv_2 >= 0, indiv_6 >= 0, indiv_7 >= 0, indiv_8 >= 0, indiv_9 >= 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "At_least_2_indiv")
    colnames(A2) <- c("cdr3aa", "group")
    
    B <-  X_CL %>% filter(pool_indiv == 0, indiv_1 > 0, indiv_2 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(B) <- c("cdr3aa", "group")
    
    C <- X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 > 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(C) <- c("cdr3aa", "group")
    
    D <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_6 > 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(D) <- c("cdr3aa", "group")
    
    E <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_6 == 0, indiv_7 > 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(E) <- c("cdr3aa", "group")
    
    G <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 > 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(G) <- c("cdr3aa", "group")
    
    H <-  X_CL %>% filter(pool_indiv == 0, indiv_1 == 0, indiv_2 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H) <- c("cdr3aa", "group")
    
    I <- X_CL %>% filter(pool_indiv > 0, indiv_1 == 0, indiv_2 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "pool_private")
    colnames(I) <- c("cdr3aa", "group")
    
    J <-  X_CL %>% filter(pool_indiv > 0, indiv_1 > 0, indiv_2 > 0, indiv_6 > 0, indiv_7 > 0, indiv_8 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "Shared_by_all")
    colnames(J) <- c("cdr3aa", "group")
    
    Indiv <- rbind(B, C, D, E, G, H)
    A2 = A2 %>% filter(!cdr3aa %in% Indiv$cdr3aa)
    A2 = A2 %>% filter(!cdr3aa %in% A$cdr3aa)
    
    CDR3aa_group = rbind(A, A2, Indiv, I, J)
    colnames(CDR3aa_group) <- c("cdr3aa", "group")
    return(CDR3aa_group)
  }
  
  if (var_cell_subset == "Teff") {
    
    A <- X_CL %>% filter(indiv_1 > 0, indiv_6 > 0, indiv_7 > 0, indiv_8 > 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "all_indiv")
    colnames(A) <- c("cdr3aa", "group")
    
    A2 <- X_CL %>% filter(indiv_1 >= 0, indiv_6 >= 0, indiv_7 >= 0, indiv_8 >= 0, indiv_9 >= 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "At_least_2_indiv")
    colnames(A2) <- c("cdr3aa", "group")
    
    B <-  X_CL %>% filter(indiv_1 > 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(B) <- c("cdr3aa", "group")
    
    D <-  X_CL %>% filter(indiv_1 == 0, indiv_6 > 0, indiv_7 == 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(D) <- c("cdr3aa", "group")
    
    E <-  X_CL %>% filter(indiv_1 == 0, indiv_6 == 0, indiv_7 > 0, indiv_8 == 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(E) <- c("cdr3aa", "group")
    
    G <-  X_CL %>% filter(indiv_1 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 > 0, indiv_9 == 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(G) <- c("cdr3aa", "group")
    
    H <-  X_CL %>% filter(indiv_1 == 0, indiv_6 == 0, indiv_7 == 0, indiv_8 == 0, indiv_9 > 0) %>%
      select(clonotypes) %>% as.data.frame() %>% mutate(group = "indiv_private")
    colnames(H) <- c("cdr3aa", "group")
    
    Indiv <- rbind(B, D, E, G, H)
    A2 = A2 %>% filter(!cdr3aa %in% Indiv$cdr3aa)
    A2 = A2 %>% filter(!cdr3aa %in% A$cdr3aa)
    
    CDR3aa_group = rbind(A, A2, Indiv)
    colnames(CDR3aa_group) <- c("cdr3aa", "group")
    return(CDR3aa_group)
  }
  
}



find_pairs <- function(x, y) {
  res <- stringdist::stringdistmatrix(x, y,
                                      method = "lv",
                                      useNames = "strings",
                                      nthread = CORES) %>%
    melt %>%
    filter(value < 2) %>%
    dplyr::select(-value)
  colnames(res) <- c("from.cdr3", "to.cdr3")
  res
}

quantile_breaks <- function(xs, n = 25) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}



fetch_stats<- function(x){
  # name=file_path_sans_ext(basename(x)) %>% str_split(pattern="/")
  print(x)
  #
  report_file<- list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN/",x),  # directory to search within
    pattern = "Stat_Clontech", 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  )
  
  if ( tools::file_ext(report_file)=="tsv"){
    all <-  read.csv(report_file, sep="\t")
  } else {
    
    all <-  read.csv(report_file)
    
  }
  print("done")
  return(all)
}

fetch_paths<- function(x){
  # name=file_path_sans_ext(basename(x)) %>% str_split(pattern="/")
  print(x)
  #
  fastq_file = list.files(
    path =paste0(choose_path,"/RS_Data/LIGAN"),  # directory to search within
    pattern = x, 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  )
  #
  fastqc_file = intersect(list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN"),  # directory to search within
    pattern = x, 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  ), 
  list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN"),  # directory to search within
    pattern = "zip", 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  ))
  #
  report_file<- intersect(list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN"),  # directory to search within
    pattern = x, 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  ), 
  list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN"),  # directory to search within
    pattern = "report.txt", 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  ))
  
  
  all <- data.frame(sample_id=x,
                    R1 =fastq_file[which(grepl( "R1",fastq_file ))],
                    R2=fastq_file[which(grepl( "R2",fastq_file ))],
                    fastqc_R1 = fastqc_file[which(grepl( "R1",fastqc_file))],
                    fastqc_R2 =fastqc_file[which(grepl( "R2",fastqc_file))],
                    report_files = report_file )
  
  
  print("done")
  return(all)
}


fetch_qc<- function(x){
  # name=file_path_sans_ext(basename(x)) %>% str_split(pattern="/")
  print(x)
  #
  report_file<- list.files(
    path =paste0(choose_path,"/RS_Processed/LIGAN/",x),  # directory to search within
    pattern = "summary_fastqc", 
    recursive = TRUE,          # search subdirectories
    full.names = TRUE          # return the full path
  )
  
  if(x=="LIGAN010"){
    
    all <-  read.table(report_file, sep=" ")
    
  } else {
    all <-  read.csv(report_file, sep="\t")
    
  }
  
  print("done")
  return(all)
}

fetch_stats_chain<- function(y){
  tr=read.csv(y, sep=",")
  chain=tools::file_path_sans_ext(basename(y)) %>% str_split(., "_")  %>% sapply(., "[", 1)
  tr_clono<- tr %>%
    filter(lib %in% freezing_fb$samplebis) %>%
    group_by(lib) %>%
    summarize(!!paste0("clonotype_number_TR",chain):=n())
  
  tr_seq<- tr %>%
    filter(lib %in% freezing_fb$samplebis) %>%
    group_by(lib) %>%
    summarize(!!paste0("sequence_number_TR",chain):=sum(count))
  
  tr_all<- merge(tr_clono, tr_seq, by="lib")
  return(tr_all)
  
}



fetch_from_freezing<- function(x){
  # Read Freezing
  freezing<-read_excel(x, sheet = "RepSeq")
  
  #keep specific columns
  freezing<- freezing %>% 
    select(A_ID_Std,'Exp. Group','Cell Sample',Lib_Exp_ID,Lib_Quant_Method,'Target genes',Lib_Labnote,Sequencer,Run_ID,
           'Sequencing direction','Read length','TCR identification method',path_processing,path_processing_full,path_filtered_full )
  
  #select samples from Batch 1 onwards
  freezing<- freezing %>% filter(`Target genes`=="TRA&TRB")
  
  #fetch [lib] and barcodes information
  aliquot=read_excel(x, sheet = "Aliquot")
  
  barcodes<- aliquot %>% 
    select(A_ID_Std,"A_Cell#","R_Quant (ng/µl)","RNAQT_used (ng)","Library_Quant (ng/ul)","Library Quantity (ng/microlitre)",'TCR Primer 2 Forward Index','TCR Primer 2 Reverse Index','Primer2 FD','Primer2 Rv',Rseq_Comment) %>%
    filter(A_ID_Std %in% freezing$A_ID_Std) %>%
    mutate(lib_conc=ifelse(!is.na(`Library_Quant (ng/ul)`) & !is.na(`Library Quantity (ng/microlitre)`), 
                           `Library Quantity (ng/microlitre)`, 
                           ifelse(is.na(`Library_Quant (ng/ul)`) & !is.na(`Library Quantity (ng/microlitre)`),
                                  `Library Quantity (ng/microlitre)`,
                                  `Library_Quant (ng/ul)`))) #%>%
  #select(-"Library_Quant (ng/ul)",-"Library Quantity (ng/microlitre)")
  
  #detect samples sequenced multiple times in different batches and duplicate rows for these samples
  # dupl<-barcodes[grepl("[(]", barcodes$`TCR Primer 2 Forward Index`),]
  dupl<-barcodes[grepl("[+]", barcodes$Rseq_Comment),]
  
  dupl<- dupl %>%
    mutate(`TCR Primer 2 Forward Index`=`TCR Primer 2 Forward Index` %>% str_split(pattern="[()]")  %>%  sapply(., "[", 2)) %>%
    mutate(`TCR Primer 2 Reverse Index`=`TCR Primer 2 Reverse Index` %>% str_split(pattern="[()]")  %>%  sapply(., "[", 2)) %>%
    mutate(Rseq_Comment=Rseq_Comment %>% str_split(pattern="[+]")  %>%  sapply(., "[", 2)) %>%
    select(-"Library_Quant (ng/ul)",-"Library Quantity (ng/microlitre)")
  
  barcodes<- barcodes %>% 
    mutate(`TCR Primer 2 Forward Index`=`TCR Primer 2 Forward Index` %>% str_split(pattern=" ")  %>%  sapply(., "[", 1)) %>%
    mutate(`TCR Primer 2 Reverse Index`=`TCR Primer 2 Reverse Index` %>% str_split(pattern=" ")  %>%  sapply(., "[", 1)) %>%
    mutate(Rseq_Comment=Rseq_Comment %>% str_split(pattern="[+]")  %>%  sapply(., "[", 1)) %>%
    mutate(lib_conc=ifelse(A_ID_Std %in% dupl$A_ID_Std , 
                           `Library_Quant (ng/ul)`,lib_conc )) %>%
    select(-"Library_Quant (ng/ul)",-"Library Quantity (ng/microlitre)")
  
  barcodes=rbind(barcodes,dupl)
  
  #add correspondance of sequencing run id
  seq_id=read_excel(x, sheet = "Correspondance_batch")
  barcodes_f<- merge(barcodes, seq_id,by.x="Rseq_Comment", by.y="Batch")
  
  #merge all primer columns
  barcodes_f<- barcodes_f %>%
    replace(is.na(.),"") %>%
    mutate(F=paste0(`TCR Primer 2 Forward Index`,`Primer2 FD`)) %>%
    mutate(R=paste0(`TCR Primer 2 Reverse Index`,`Primer2 Rv`)) 
  
  
  #                 
  corrbarcodes=read_excel(x, sheet = "5Racebarcodes")
  
  corrbarcodes=left_join(as.data.frame(barcodes_f),as.data.frame(corrbarcodes),
                         by=c("F"="Forward", "R"="Reverse"))
  
  # corrbarcodes_1<- corrbarcodes_1 %>% select(A_ID_Std,barcode,sequencing_run) %>%
  #                   mutate(barcode=as.character(barcode)) %>% 
  #                   replace(is.na(.),0) %>%
  #                   filter(barcode!=0)
  # 
  # 
  # corrbarcodes_2=left_join(as.data.frame(barcodes_f),as.data.frame(corrbarcodes), by=c("Primer2 FD"="Forward","Primer2 Rv"="Reverse"))
  # corrbarcodes_2<- corrbarcodes_2 %>% select(A_ID_Std,barcode,sequencing_run) %>%
  #                   mutate(barcode=as.character(barcode)) %>% 
  #                   replace(is.na(.),0) %>%
  #                   filter(barcode!=0)
  # 
  # 
  # both_barcodes<- rbind(corrbarcodes_1,corrbarcodes_2)
  # both_barcodes$A_ID_Std[which(duplicated(both_barcodes$A_ID_Std))]
  # barcodes$A_ID_Std[which(!barcodes$A_ID_Std %in% both_barcodes$A_ID_Std)]
  
  # add a step to merge freezing and barcdes_f and then use it for the next step
  freezing_f<- merge(freezing,corrbarcodes,by.x=c("A_ID_Std","Run_ID"), by.y=c("A_ID_Std","sequencing_run"))                     
  
  return(freezing_f)
}

