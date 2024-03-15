source("./Functions.R")
CORES = 4

#Figure S1 ####

all_data<- all_data %>% mutate(cell_subset=factor(cell_subset, levels=c("CD8","Teff", "nTregs","amTregs")))
#a
ggplot(all_data, aes(x =cell_number, y =RNA))+
  geom_point(aes(fill=`Per base sequence quality`, shape=organ), color="black")+
  facet_grid(~ cell_subset, scales="free")+
  smplot2::sm_statCorr(color = 'black', corr_method = 'spearman', linetype = 'dashed',text_size = 3)+
  scale_shape_manual(values=c(21,22,24))+
  theme_Publication()+
  ggplot2::theme( axis.text = element_text(size=7))+
  scale_x_continuous(labels = function(x) format(x, scientific =TRUE))+
  labs(title = "Cell number vs. [RNA]")

#b
ggplot(all_data_m[grepl("clonotype_", all_data_m$variable) ,], 
       aes(x=organ, y=as.numeric(value)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(fill=`Per base sequence quality`),width=0.2,shape=21)+
  facet_grid(variable~ cell_subset)+
  theme_Publication()+
  ggplot2::theme( axis.text = element_text(size=7))+
  ylab('total number of clonotypes')+
  xlab("")+
  labs(title = "Total number of clonotypes")

#c
correlation_stats=all_data %>%
  filter(`Per base sequence quality`=="PASS") %>%
  reshape2::melt(id=c("sample_id","cell_subset","organ","RNA","RNAQT_used (ng)", "cell_number","cDNA")) %>%
  filter(grepl(paste(c("clonotype_number_","sequence_number_"), collapse="|"), variable))

correlation_stats$value<- as.numeric(correlation_stats$value) 


sp <- ggscatter(correlation_stats %>%  filter(grepl("clonotype_number_", variable)), x = "RNA", y = "value",
                fill="cell_subset",
                add = "reg.line",  # Add regressin line
                add.params = list( color="cell_subset"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                ylab="# of clonotypes",
                xlab="[RNA]",
                facet.by = "variable",
                palette = "jco",
                shape=21
                
)
# Add correlation coefficient
sp + stat_cor(method = "spearman", aes(color=cell_subset))+
  theme_classic()+theme(plot.margin = unit(c(1,2,1,2) ,"cm"))


#d
sp <- ggscatter(correlation_stats %>%  filter(grepl("clonotype_number_", variable)), x = "cell_number", y = "value",
                fill="cell_subset",
                add = "reg.line",  # Add regressin line
                add.params = list(color="cell_subset"), # Customize reg. line
                conf.int = FALSE, # Add confidence interval
                ylab="# of clonotypes",
                xlab="cell_number",
                facet.by = "variable",
                palette = "jco",
                shape=21
                
)
# Add correlation coefficient
sp + stat_cor(method = "spearman", aes(color=cell_subset))+
  theme_classic()+theme(plot.margin = unit(c(1,2,1,2) ,"cm"))

#e
indiv_multiple$quantity<- indiv_multiple$Volume*indiv_multiple$R_Quant..ng.µl.
ggplot(indiv_multiple, aes(x =cell_number, y =quantity))+
  geom_point(aes(fill=as.factor(poolplus)),shape=21, size=2.5,color="black")+
  facet_grid(~ cell_subset, scales="free")+
  smplot2::sm_statCorr(color = 'black', corr_method = 'spearman', linetype = 'dashed',text_size = 3)+
  theme_Publication()+
  ggplot2::theme( axis.text = element_text(size=7))+
  scale_size_area(max_size=4)+
  scale_fill_brewer(type="seq", palette="YlOrRd")+
  ylab("Total RNA quantity")

#f
df_Rarefaction_res<- fortify(Rarefaction_res)
df_Rarefaction_res <- separate(df_Rarefaction_res, Assemblage, into=c("poolplus", "sample_id", "cell_subset", "chain"), sep="%" )

df.point <- df_Rarefaction_res[which(df_Rarefaction_res$Method=="Observed"),]
df.line <- df_Rarefaction_res[which(df_Rarefaction_res$Method !="Observed"),]

df_Rarefaction_res$poolplus=factor(df_Rarefaction_res$poolplus, levels=c("1", "2","3","6" ,"7", "8" ,"9" , "10"))
df_Rarefaction_res$cell_subset=factor(df_Rarefaction_res$cell_subset,levels=c("CD8","Teff","nTregs","amTregs"))

df.point$poolplus=factor(df.point$poolplus, levels=c("1", "2","3","6" ,"7", "8" ,"9" , "10"))
df.point$cell_subset=factor(df.point$cell_subset,levels=c("CD8","Teff","nTregs","amTregs"))

df.line$poolplus=factor(df.line$poolplus, levels=c("1", "2","3","6" ,"7", "8" ,"9" , "10"))
df.line$cell_subset=factor(df.line$cell_subset,levels=c("CD8","Teff","nTregs","amTregs"))

ggplot(df_Rarefaction_res, aes(x=x, y=y)) +
  geom_line(data = df.line %>% filter(Method == "Rarefaction"), aes(x=x, y=y, colour = factor(poolplus), group = sample_id), linetype = "solid", lwd=0.8,  alpha = 0.8) +
  geom_line(data = df.line %>% filter(Method == "Extrapolation"), aes(x=x, y=y, colour = factor(poolplus), group = sample_id), linetype = "dashed", lwd=0.8, alpha = 0.8) +
  geom_point(data=df.point, aes(fill = factor(poolplus)), colour = "black", shape = 21, alpha = 0.9) +
  labs(x="Sample size", y="Number of unique clonotype") +
  facet_grid(chain~cell_subset, scales = "free")+
  scale_fill_brewer(type="seq", palette="YlOrRd")+
  scale_color_brewer(type="seq", palette="YlOrRd")+
  theme_Publication()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


#Figure S2####
# Create the summary table
summary <- dt_filtered %>%
            dplyr::group_by(filename, Mice_ID, chain, cell_subset, Cell_Sample, cell_number) %>%
            dplyr::mutate(Nbr_clonotypes = n_distinct(clonotypes)) %>%
            dplyr::mutate(Nbr_sequences = sum(cloneCount)) %>%
            dplyr::mutate(Nbre_seq_per_cell = Nbr_sequences/cell_number) %>%
            dplyr::mutate(Nbre_clo_per_sequences = Nbr_clonotypes/Nbr_sequences) %>%
            dplyr::distinct(filename, cell_subset, chain, Nbr_clonotypes, Nbr_sequences, Nbre_seq_per_cell,Nbre_clo_per_sequences, Cell_Sample, cell_number)

##colors
nejm_color <- c("#e18727",  "#20854e","#bc3c29", "#7876b1", "#0072b5")


##a
fout <- c("tripod-64-1923_R1.txt","tripod-64-1927_R1.txt","tripod-64-1931_R1.txt",
          "tripod-64-1917_R1.txt","tripod-64-1921_R1.txt","tripod-64-1925_R1.txt",
          "tripod-64-1929_R1.txt","tripod-64-1932_R1.txt","tripod-93-2846_R1.txt")

summary_temp<- summary %>% 
  filter(!filename %in% fout ) %>%
  select(-Mice_ID,-cell_number,-Nbre_seq_per_cell,-Nbre_clo_per_sequences) %>%
  reshape2::melt(by=c("filename","chain","cell_subset","Cell_Sample"))

summary_temp$cell_subset = if_else(summary_temp$cell_subset == "CD4+", "Teff", summary_temp$cell_subset)

summary_temp$cell_subset <- factor(summary_temp$cell_subset, levels = c("CD8+", "Teff", "nTreg", "amTreg"))
summary_temp$Cell_Sample <- factor(summary_temp$Cell_Sample, levels = c("indiv", "pool"))

stat.test <- summary_temp %>% 
  filter(variable == "Nbr_sequences") %>%
  group_by(cell_subset, chain) %>%
  rstatix::wilcox_test(value ~ Cell_Sample) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>% 
  rstatix::add_xy_position()

## Graph
ggpubr::ggboxplot(data = summary_temp %>% filter(variable == "Nbr_sequences"),
                  x = "Cell_Sample", y = "value", alpha=.5,  fill = "Cell_Sample", color = "Cell_Sample", outlier.shape = NA,
                  facet.by = c("chain","cell_subset"), scales ="free_y")+
  geom_jitter( aes( fill=Cell_Sample),size=1.5,color = "black", shape=21,alpha=.4, position = position_jitterdodge(.2))+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  scale_y_continuous( "Number of sequences", expand = expansion(mult = c(0.1, 0.1))) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01, 
    hide.ns = FALSE)

##b
stat.test <- summary_temp %>% 
  filter(variable == "Nbr_clonotypes") %>%
  group_by(cell_subset, chain) %>%
  rstatix::wilcox_test(value ~ Cell_Sample) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>% 
  rstatix::add_xy_position()

## Graph
ggpubr::ggboxplot(data = summary_temp %>% filter(variable == "Nbr_clonotypes"),
                  x = "Cell_Sample", y = "value", alpha=.5,  fill = "Cell_Sample", color = "Cell_Sample", outlier.shape = NA,
                  facet.by = c("chain","cell_subset"), scales ="free_y")+
  geom_jitter( aes( fill=Cell_Sample),size=1.5,color = "black", shape=21,alpha=.4, position = position_jitterdodge(.2))+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  scale_y_continuous( "Number of clonotypes", expand = expansion(mult = c(0.1, 0.1))) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01, 
    hide.ns = FALSE)



#Figure 1####
##a
n_TR = max(summary$Nbr_sequences)

temp <- dt_filtered %>% named_group_split(Cell_Sample, cell_subset, chain, sep="")
temp <- mclapply(dt_filtered, data_modelling, row = "group2", col = "clonotypes", value = "cloneCount", mc.cores = 4)
temp <- data_modelling(dt_filtered, row = "group2", col = "clonotypes", value = "cloneCount")


Rarefaction_res <- iNEXT(temp, q=0, datatype="abundance", endpoint = n_TR)
Rarefaction_res = rbindlist(lapply(Rarefaction_res, setDT, keep.rownames = TRUE))
mod <- setnames(mod, old = "rn", new = "group2")

Rarefaction_res <- ggplot2::fortify(Rarefaction_res,1)
Rarefaction_res$A = Rarefaction_res$site
Rarefaction_res <- separate(Rarefaction_res, A, into=c("Cell_Sample", "Mice_ID", "cell_subset", "chain", "Clin_Group"), sep="%" )

Rarefaction_res$cell_subset = if_else(Rarefaction_res$cell_subset == "CD8", "CD8+",
                                      if_else(Rarefaction_res$cell_subset == "Teff", "Teff",
                                              if_else(Rarefaction_res$cell_subset == "amTregs", "amTreg",
                                                      if_else(Rarefaction_res$cell_subset == "nTregs", "nTreg","x"))))

df.point <- Rarefaction_res[which(Rarefaction_res$method=="observed"),]
df.line <- Rarefaction_res[which(Rarefaction_res$method!="observed"),]
df.line$method <- factor(df.line$method,
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))


Rarefaction_res$Cell_Sample <- factor(Rarefaction_res$Cell_Sample, levels = c("indiv", 'pool'))
Rarefaction_res$cell_subset <- factor(Rarefaction_res$cell_subset, levels = c("CD8+", 'Teff', 'nTreg', "amTreg"))
df.line$Cell_Sample <- factor(df.line$Cell_Sample, levels = c("indiv", 'pool'))
df.point$Cell_Sample <- factor(df.point$Cell_Sample, levels = c("indiv", 'pool'))
df.point$cell_subset <- factor(df.point$cell_subset, levels = c("CD8+", 'Teff', 'nTreg', "amTreg"))
df.line$cell_subset <- factor(df.line$cell_subset, levels = c("CD8+", 'Teff', 'nTreg', "amTreg"))

## Graph
ggplot(Rarefaction_res, aes(x=x, y=y)) +
  geom_line(data = df.line %>% filter(method == "interpolation"), aes(x=x, y=y, colour = Cell_Sample, group = Mice_ID), linetype = "solid", lwd=0.8, show.legend = F, alpha = 0.8) +
  geom_line(data = df.line %>% filter(method == "extrapolation"), aes(x=x, y=y, colour = Cell_Sample, group = Mice_ID), linetype = "dashed", lwd=0.8, show.legend = F, alpha = 0.8) +
  geom_point(data=df.point, aes(fill = Cell_Sample), colour = "black", shape = 21, alpha = 0.9, show.legend = T) +
  labs(x="Sample size", y="Number of unique clonotype") +
  scale_y_continuous(labels = function(x) format(x, scientific = T))+
  facet_grid(chain~cell_subset, scales = "free")+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

##S2 c
stat_rarefaction = Rarefaction_res %>%
  group_by(cell_subset, chain, Cell_Sample, Mice_ID) %>%
  summarise(max_y = max(y))

stat_rarefaction <- stat_rarefaction %>%
  group_by(cell_subset, chain) %>%
  rstatix::wilcox_test(max_y ~ Cell_Sample) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj")

stat_rarefaction <- Rarefaction_res %>%
  group_by(cell_subset, chain, Cell_Sample, Mice_ID) %>%
  mutate(AUC = MESS::auc(x, y, type = 'spline')) %>%
  select(Cell_Sample, AUC, chain, cell_subset) %>% distinct()


## Graph
stat_rarefaction %>% select(Cell_Sample, AUC, chain, cell_subset) %>% distinct() %>%
  ggplot(
    aes(Cell_Sample, AUC, fill = Cell_Sample))+
  geom_jitter(color = "black", shape =20, size = 3, width = 0.15, height = 0, alpha = 0.5, show.legend = F) +
  geom_boxplot(aes(color = Cell_Sample), fill = NA, alpha = 1, outlier.shape = NA) +
  facet_grid(chain ~ cell_subset)+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  ylab("Rarefaction (AUC)") +
  ggpubr::stat_compare_means(label= "p.signif", comparisons = my_comparisons, p.adjust.methods = "Holm") +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) + 
  theme_Publication() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

##b
temp <- dt_filtered %>% 
       named_group_split(Cell_Sample, cell_subset, chain, sep="%")

temp <- mclapply(temp, data_modelling, row = "group2", col = "clonotypes", value = "cloneCount", mc.cores = 3)

mod <- mclapply(temp, vegan::renyi, mc.cores = 4)
mod = rbindlist(lapply(mod, setDT, keep.rownames = TRUE))

setnames(mod, old = "rn", new = "group2")

summary_temp$cell_subset<- gsub("Teff", replacement = "CD4+",summary_temp$cell_subset)
summary_temp$group2 = paste(summary_temp$Cell_Sample, summary_temp$Mice_ID, summary_temp$cell_subset, summary_temp$chain, "B6_Yg", sep = "%")
summary_temp<- summary_temp %>% select(-variable, -value) %>% distinct()

mod <- merge(mod, summary_temp[,c( "Mice_ID","cell_subset", "Cell_Sample", "chain", "group2")], by = "group2")
mod <- melt(mod)
mod2 <- Rmisc::summarySE(mod,  measurevar = "value", 
                         groupvars=c("variable", "cell_subset","chain", "Cell_Sample"))
rm(temp)
mod2$Cell_Sample <- factor(mod2$Cell_Sample, levels = c("indiv", 'pool'))
mod$Cell_Sample <- factor(mod$Cell_Sample, levels = c("indiv", 'pool'))
mod2$cell_subset <- factor(mod2$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))
mod$cell_subset <- factor(mod$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

## Graph 1
ggplot(mod2,
       aes(x=variable, y=value, group=Cell_Sample, color=Cell_Sample, fill=Cell_Sample)) +
  geom_ribbon(aes(ymin=value-se, ymax=value+se), alpha=0.2, color=NA)+
  geom_line(alpha = 0.7) +
  geom_point(shape=24,size=2, color="black", alpha = 0.7)+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  ylab("Renyi entropy") +
  xlab("alpha")+
  facet_grid(chain~cell_subset)+
  scale_y_continuous(expand = expand_scale(mult = c(0, .15))) + 
  theme_Publication()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

stat.test <- data.frame()
## Graph 2
for(v in unique(mod$variable)) {
  
  temp1 <- mod %>% filter(variable == v) %>%
    group_by(variable, cell_subset, chain) %>%
    rstatix::wilcox_test(value ~ Cell_Sample) %>%
    rstatix::adjust_pvalue(method = "holm") %>%
    rstatix::add_significance("p.adj")
  
  stat.test <- rbind(stat.test, temp1)
}


stat.test <- stat.test %>% filter(variable %in% c(0, 1, 8, Inf))

##S2 d
X_CL <- dt_filtered

X_CL$freq_bin<-cut(X_CL$freq_clo, c(0.000001, 0.00001, 0.0001,0.001,0.01, 0.1,1, 100), include.lowest=T)

X_CL <- X_CL %>%
  group_by(Mice_ID, chain, cell_subset, Clin_Group, Cell_Sample, freq_bin) %>%
  #summarise(count = n()) %>%
  summarise(count = sum(freq_clo)) %>%
  ungroup() %>%
  group_by(Mice_ID, chain, cell_subset, Clin_Group, Cell_Sample) %>%
  mutate(percentage = count/sum(count))


X_CL$Cell_Sample <- factor(X_CL$Cell_Sample, levels = c("indiv", 'pool'))
X_CL$cell_subset <- factor(X_CL$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

X_occ_2 <- Rmisc::summarySE(X_CL, measurevar = "count", groupvars = c("Clin_Group", "cell_subset", "chain", "Cell_Sample", "freq_bin"))
X_occ_2$Cell_Sample <- factor(X_occ_2$Cell_Sample, levels = c("indiv", 'pool'))
X_occ_2$cell_subset <- factor(X_occ_2$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

ggplot2::ggplot(X_occ_2, aes(x = freq_bin, y = ifelse(is.na(count), 0, count), fill = Cell_Sample)) +
  facet_grid(chain ~ cell_subset, scales = "fixed", space = "fixed")+
  geom_bar(stat="identity", color = "black", width=0.8, position = position_dodge(preserve = "single"), show.legend = T) +
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,
                position=position_dodge(.9))+
  #geom_density(stat = "identity", alpha = 0.3, aes(group = Cell_Sample, color = Cell_Sample), fill = NA)+
  scale_color_manual(values = nejm_color)+
  scale_fill_manual(values = nejm_color)+
  scale_y_continuous("Cumulative frequency", expand = expansion(mult = c(0, .15))) + 
  theme_Publication() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


##d
X_CL$freq_bin <- as.factor(X_CL$freq_bin)

ggplot(X_CL, aes(fill=as.factor(freq_bin), y=percentage, x=Mice_ID)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(chain~cell_subset+Cell_Sample, scales = "free") +
  labs(x = "Samples", y = "Cumulative frequency") +
  scale_fill_brewer(palette="RdYlBu", direction=-1)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 6))


temp50 <- X_CL %>% group_by(Mice_ID, chain, cell_subset, Cell_Sample) %>%
  complete(freq_bin, fill = list(count = 0, percentage = 0))

stat_res = data.frame()
for(fb in unique(temp50$freq_bin)) {
  temp1 <- temp50 %>% filter(freq_bin == fb)
  temp1 <- compare_means(count~Cell_Sample, data = temp1,
                         group.by = c("freq_bin", "chain", "cell_subset")) %>%
    rstatix::adjust_pvalue(method = "holm") %>%
    rstatix::add_significance("p.adj")
  stat_res <- rbind(stat_res, temp1)
}


#Figure 2 ####
##a

X_CL = dt_filtered

X_CL <- X_CL %>%
  group_by(filename, chain) %>%
  arrange(desc(cloneCount)) %>%
  mutate(clone_ID = row_number()) %>%
  mutate(quartile = cut(clone_ID, breaks = 4, labels = paste0("Q", 1:4)))

summary_quartile <- X_CL %>%
  group_by(filename, chain, Cell_Sample, Clin_Group, cell_subset, quartile) %>%
  summarise(n_distinct(clonotypes))
# 
# X_CL_Q1 <- X_CL %>% filter(quartile == "Q1")
# X_CL_Q1 <- X_CL_Q1 %>% named_group_split(chain,Cell_Sample, Clin_Group, cell_subset, quartile, sep = "_")
# X_CL_Q1 <- lapply(X_CL_Q1, dkast, gene = "clonotypes", value = "cloneCount")

X_CL2 <- X_CL %>% named_group_split(chain, Clin_Group, cell_subset, Cell_Sample, sep="%")
X_CL2 <- lapply(X_CL2, dkast, gene = "clonotypes", value = "cloneCount")

MH_res3 <- mclapply(X_CL2, Compute_MH, method = "ji", out="df", mc.cores = 4) ## method mh or ji

MH_res2 <- rbindlist(MH_res3, idcol = "group")
MH_res2 <- separate(MH_res2, Var1, into = c("Cell_Sample1", "Mice_ID1", "cell_subset1", "chain1", "Clin_Group1"), sep = "%")
MH_res2 <- separate(MH_res2, Var2, into = c("Cell_Sample2", "Mice_ID2", "cell_subset2", "chain2", "Clin_Group2"), sep = "%")
MH_res2 <- separate(MH_res2, group, into = c("V1", "V2", "V3", "Quantile"), sep = "%")
MH_res2 <- MH_res2 %>% select(-V1, -V2, -V3)

MH_res2 = MH_res2 %>%
  filter(cell_subset1 == cell_subset2) %>%
  filter(chain1 == chain2) %>%
  filter(Clin_Group1 == Clin_Group2) %>%
  filter(Cell_Sample1 == Cell_Sample2) 

MH_res2 <- MH_res2[!duplicated(data.frame(t(apply(MH_res2,1,sort)))),]
MH_res2$cell_subset1 <- factor(MH_res2$cell_subset1, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

## Graph 1
stat.test <- MH_res2 %>%
  group_by(chain1, cell_subset1) %>%
  rstatix::wilcox_test(value ~ Cell_Sample1) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>% rstatix::add_xy_position()


ggpubr::ggboxplot(data = MH_res2,
                  x = "Cell_Sample1", y = "value", alpha=.5,  fill = "Cell_Sample1", color = "Cell_Sample1", outlier.shape = NA,
                  facet.by = c("chain1","cell_subset1"), scales ="free_y")+
  geom_jitter( aes( fill=Cell_Sample1),size=1.5,color = "black", shape=21,alpha=.4, position = position_jitterdodge(.2))+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), "Morisita Horn similarity score", expand = expansion(mult = c(0.1, 0.1))) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01, 
    hide.ns = TRUE)

##b
library(GGally)
X_CL = dt_filtered


X_CL <- X_CL %>%
  group_by(filename, chain) %>%
  arrange(desc(cloneCount)) %>%
  mutate(clone_ID = row_number()) %>%
  mutate(quartile = cut(clone_ID, breaks = 4, labels = paste0("Q", 1:4)))

X_CL_Q1 <- X_CL %>% filter(quartile == "Q1")
X_CL_Q1 <- X_CL_Q1 %>% named_group_split(chain, sep="%")
X_CL_Q1 <- lapply(X_CL_Q1, dkast, gene = "clonotypes", value = "freq_clo")

data_cor <- lapply(X_CL_Q1, function(x) x <- cor(x))

files_CL <- data_cor$`TRB` # or TRA
column_group <- X_CL %>% ungroup() %>% distinct(group2, cell_subset, Cell_Sample) %>% as.data.frame()
rownames(column_group) <- column_group$group2
column_group <- as.data.frame(column_group)
column_group <- column_group %>% select(-group2)

### Heatmap preparation
mat_colors <- list(
  cell_subset = c("CD8+" = "#1F78B4", "CD4+" = "#bc3c29", "nTreg" = "#72BF5B", "amTreg" = "#7876b1"),
  Cell_Sample = c(pool = "#20854e", indiv = "#e18727"))

#nejm_color <- c("#e18727",  "#20854e","#bc3c29", "#7876b1", "#0072b5")

### Colours choice
break_ref <- colnames(files_CL)[max.col(files_CL, ties.method="first")]
mat_breaks <- quantile_breaks(files_CL[,break_ref[1]], n = 2000)

str(files_CL)
str(column_group)
setdiff(rownames(column_group), colnames(files_CL))

column_group <- column_group %>% filter(rownames(column_group) %in% colnames(files_CL))

### Plot the heatmap
pheatmap::pheatmap(files_CL,
                   #breaks = mat_breaks,
                   #color = viridis::turbo(length(mat_breaks) - 1), 
                   color = viridis::cividis(n = 50), 
                   annotation_colors = mat_colors,
                   
                   #cluster_cols=as.hclust(col_dend),
                   show_rownames = F, show_colnames = F,
                   #scale = "column",
                   border_color = NA,
                   # cutree_cols = 4, cutree_rows = 4,
                   main = "TRB", 
                   cluster_rows = T, #cluster_cols = T,
                   clustering_distance_rows="correlation", clustering_distance_cols ="correlation",
                   clustering_method="ward.D2",
                   annotation_names_col = T,
                   annotation_row = column_group, annotation_col = column_group)

##c
plot <- X_CL_Q1 %>%
  group_by(Cell_Sample,cell_subset) %>%
  mutate(id2=ifelse(grepl("_1",Mice_ID), "S1", "S2")) %>%
  filter(Mice_ID == "indiv_1" | Mice_ID == "indiv_9" | Mice_ID == "pool_1" | Mice_ID == "pool_2"  )
plot_d<- dcast(plot, clonotypes+chain+cell_subset+Cell_Sample~id2, value.var="cloneCount")
plot_d[is.na(plot_d)]<-0
plot_d$mean<- (plot_d$S1+plot_d$S2)/2
plot_d$cell_subset <- factor(plot_d$cell_subset, levels = c("CD8+", "CD4+", "nTreg", "amTreg"))


p_ <- GGally::print_if_interactive
p_(ggally_points(plot_d %>% filter(chain=="TRB") , # or TRA
                 mapping = ggplot2::aes_string(x = "S1", y = "S2",  
                                               fill="log2(mean)"),#shape="Cell_Sample",
                 alpha = 0.7, shape=21, size=3))+
  facet_grid(.~cell_subset+Cell_Sample)+
  scale_x_log10() +
  scale_y_log10()+
  theme_Publication()+
  # scale_shape_manual(values=c(21,21))+
  scale_fill_viridis(limits = c(0.5, 12.1))+ 
  # theme(text = element_text(family = "Arial"))
  theme(legend.position = "none",text = element_text(family = "Arial"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

##d
X_CL_Q1 <- X_CL %>% filter(quartile == "Q1")

X_CL_Q1 = X_CL_Q1 %>% 
  group_by(Sample_ID_Std,chain,cell_subset,Cell_Sample,length_cdr3aa) %>%
  summarise(count=n()) %>%
  group_by(Sample_ID_Std,chain,cell_subset,Cell_Sample) %>%
  dplyr::mutate(percent=count/sum(count))

data_filtered_summ<-Rmisc::summarySE(X_CL_Q1, measurevar="percent",
                                     groupvars=c("chain", "cell_subset","Cell_Sample","length_cdr3aa"),conf.interval = 0.95)

X_CL_Q1$cell_subset <-  factor(X_CL_Q1$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))
data_filtered_summ$cell_subset <-  factor(data_filtered_summ$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

ggplot(X_CL_Q1, aes(x=length_cdr3aa, y=percent, color = Cell_Sample)) +
  # geom_jitter( position = position_jitter(0.2), alpha=.5) + 
  geom_line(aes(group = Cell_Sample), data =data_filtered_summ ) +
  geom_errorbar(
    aes(ymin = percent-sd, ymax = percent+sd),
    data = data_filtered_summ , width = 0.15) +
  facet_grid(chain~cell_subset)+
  geom_point(data = data_filtered_summ , aes(fill=Cell_Sample), size = 1, shape=21)+
  scale_color_manual(values =c(pool = "#20854e", indiv = "#e18727"))+
  scale_fill_manual(values =c(pool = "#20854e", indiv = "#e18727"))+
  ggpubr::stat_compare_means(aes(group = Cell_Sample), label = "p.signif", hide.ns = TRUE, label.y = 0.33)+
  theme_Publication()

##e
cdr3aa_beta= as.data.frame(unique(X_CL_Q1$cdr3aa[X_CL_Q1$chain=="TRB"]))

cdr3aa_alpha= as.data.frame(unique(X_CL_Q1$cdr3aa[X_CL_Q1$chain=="TRA"]))

# run olga-compute_pgen --mouseTRB -i

# file_list <- list.files(path = "./files/", pattern = "pgen", include.dirs = T)
pgen_filename = file_list
names(file_list) <- file_list
# file_list <- lapply(paste0( "./files/", file_list), fread)
file_list <- mapply(cbind, file_list, "filename" = pgen_filename, SIMPLIFY = F)
pgen_dt = rbindlist(file_list)
colnames(pgen_dt)<- c("cdr3aa","pgen","filename")

X_CL_Q1 <- X_CL %>% filter(quartile == "Q1")

networks_pgen=merge(X_CL_Q1,pgen_dt[,-3], by = "cdr3aa")
med_res<- networks_pgen %>%
  filter(pgen>10e-20) %>%
  mutate(cell_subset=factor(cell_subset,levels = c("CD8+", "CD4+", "nTreg", "amTreg") )) %>%
  dplyr::group_by(chain,cell_subset,Cell_Sample) %>%
  dplyr::summarize(med=median(pgen))

networks_pgen %>%
  filter(pgen>10e-20) %>%
  mutate(cell_subset=factor(cell_subset,levels = c("CD8+", "CD4+", "nTreg", "amTreg") )) %>%
  ggplot( aes(x=pgen, color=Cell_Sample, group=Mice_ID))+
  geom_density()+
  facet_grid(chain~cell_subset)+
  scale_x_log10()+
  geom_vline(aes(xintercept=med, col=Cell_Sample), data=med_res, linetype="dashed")+
  scale_color_manual(values =c(pool = "#20854e", indiv = "#e18727"))+
  theme_Publication() 

#Figure S3####

# *********************************************** #
#### MetaTCR creation ####
# *********************************************** #

# Method 1####
# Set size equal to pool_indiv sample (by chain and cell_subset)
size_dt <- dt %>%
  dplyr::filter(Mice_ID == "pool_indiv") %>%
  dplyr::group_by(cell_subset, chain) %>%
  dplyr::summarise(Sequences_number = round(sum(cloneCount))) %>%
  dplyr::filter(cell_subset %in% c("amTreg", "nTreg")) %>%
  dplyr::mutate(Nb_for_5 = round(Sequences_number/5)) %>%
  dplyr::mutate(Nb_for_6 = round(Sequences_number/6))

Mice_size <- dt %>%
  dplyr::filter(Cell_Sample == "indiv") %>%
  dplyr::distinct(Mice_ID, cell_subset, chain, nb_sequences)

# Add this information into dt of interest
Mice_size <- merge(Mice_size, size_dt, by = c("cell_subset", "chain"))
Mice_size <- ungroup(Mice_size)

Mice_size$keep <- ifelse(Mice_size$cell_subset == "amTreg" & Mice_size$nb_sequences > Mice_size$Nb_for_5, "Yes",
                         ifelse(Mice_size$cell_subset == "nTreg" & Mice_size$nb_sequences > Mice_size$Nb_for_5, "Yes", "No"))

dt <- dt %>%
  dplyr::filter(cell_subset %in% c("amTreg", "nTreg"))

dt_indiv = dt %>%
  dplyr::filter(Cell_Sample == "indiv") %>%
  dplyr::filter(!(cell_subset == "amTreg" & Mice_ID == "indiv_6"))

dt_metaTCR = data.frame()

for (i in 1:10) {
  print(i)
  for (ch in unique(dt_indiv$chain)) {
    for (cs in unique(dt_indiv$cell_subset)) {
      for (mID in unique(dt_indiv$Mice_ID)) {
        
        if (cs == "amTreg" & mID == "indiv_6") {
          print(paste0("/ NA /"))
        }  else {      
          
          size_temp <- size_dt %>%
            dplyr::filter(chain == ch) %>%
            dplyr::filter(cell_subset == cs)
          
          if (cs == "amTreg") {
            size_temp <- size_temp %>% ungroup() %>%
              dplyr::select(Nb_for_5) %>% unique() %>% unlist()
          }
          
          if (cs == "nTreg") {
            size_temp <- size_temp %>% ungroup() %>%
              dplyr::select(Nb_for_6) %>% unique() %>% unlist()
          }
          cat(ch, cs, mID, size_temp, " / NEXT / ", sep = " / ")
          
          dt_temp = dt_indiv %>%
            dplyr::filter(chain == ch) %>%
            dplyr::filter(cell_subset == cs) %>%
            dplyr::filter(Mice_ID == mID) %>%
            dplyr::slice(rep(1:n(), cloneCount)) %>%
            dplyr::select(-cloneCount) %>%
            dplyr::mutate(cloneCount = 1) %>%
            sample_n(size = size_temp) %>% 
            dplyr::group_by_all() %>%
            dplyr::summarise(cloneCount = sum(cloneCount)) %>% ungroup()
          
          sum(dt_temp$cloneCount) == size_temp
          
          dt_temp$iteration = i
          dt_metaTCR = rbind(dt_metaTCR, dt_temp)
          
        }
      }
      gc(reset = T)
    }
  }
}

dt_metaTCRsave <- dt_metaTCR
dt_metaTCR$Mice_ID <- paste0("MetaTCR", dt_metaTCR$iteration)

dt_metaTCR_summary <- dt_metaTCR %>%
  dplyr::group_by(cell_subset, chain, Mice_ID) %>%
  dplyr::summarise(Sequences_number = round(sum(cloneCount)))

dt_metaTCR <- dt_metaTCR %>% select(-iteration)

##aggregate counts
dt_metaTCR_agg<- dt_metaTCR %>% 
  select(chain,clonotypes,Mice_ID,Clin_Group,sex,cell_subset, cloneCount) %>%
  group_by(chain,clonotypes,Mice_ID,Clin_Group,sex,cell_subset) %>%
  summarize(cloneCount=sum(cloneCount)) %>%
  group_by(chain,Mice_ID,Clin_Group,sex,cell_subset) %>%
  mutate(freq_clo=cloneCount/sum(cloneCount)*100) %>%
  mutate(filename=Mice_ID) %>%
  mutate(Cell_Sample="MetaTCR 1") %>%
  ungroup()%>%
  select(filename,cloneCount,chain,clonotypes,freq_clo,Mice_ID,cell_subset,Cell_Sample)


t1<- dt_metaTCR_agg  %>%
  group_by(cell_subset, chain, filename,Cell_Sample) %>%
  summarize(count=n()) %>%
  group_by(cell_subset, chain,Cell_Sample) %>%
  summarise( 
    n=n(),
    mean=mean(count),
    sd=sd(count)
  ) %>%
  mutate( se=sd/sqrt(n))  

# Method 3 ####
keep <- dt %>%
  filter(Cell_Sample == "indiv" & grepl("Treg",cell_subset ) ) %>%
  mutate(count=1) %>%
  reshape2::dcast(., cell_subset+chain+clonotypes~filename) %>%
  mutate(sum=rowSums(select_if(., is.numeric),na.rm = TRUE)) %>%
  filter(sum==1) %>%
  select(cell_subset, chain, clonotypes)

dt_bis <- dt %>%
  filter(Cell_Sample == "indiv" & grepl("Treg",cell_subset ) ) %>%
  merge(., keep, by=c("cell_subset", "chain", "clonotypes")) %>%
  group_by(filename, cell_subset,chain) %>%
  mutate(freq_clo=cloneCount/sum(cloneCount)*100) %>%
  mutate(Cell_Sample="MetaTCR 3") %>%
  mutate(Mice_ID=filename) %>%
  mutate(filename="MetaTCR 3") %>%
  select(filename, cloneCount,chain,clonotypes, freq_clo,Mice_ID, cell_subset, Cell_Sample) 

t3<- dt_bis %>%
  group_by(cell_subset, chain,Cell_Sample) %>%
  summarize(count=n()) %>%
  mutate( 
    n=1,
    mean=count,
    sd=NA,
    se=NA
  ) %>%
  select(-count)

# Method 2 ####
dt_bis2 <- dt %>%
  filter(Cell_Sample == "indiv" & grepl("Treg",cell_subset ) ) %>%
  group_by(cell_subset, chain,  clonotypes) %>%
  summarize(cloneCount=sum(cloneCount)) %>%
  group_by( cell_subset,chain) %>%
  mutate(freq_clo=cloneCount/sum(cloneCount)*100) %>%
  mutate(Cell_Sample="MetaTCR 2") %>%
  mutate(Mice_ID="MetaTCR 2") %>%
  mutate(filename=Mice_ID) %>%
  select(filename, cloneCount,chain,clonotypes, freq_clo, Mice_ID,cell_subset, Cell_Sample) 

t2<- dt_bis2 %>%
  group_by(cell_subset, chain,Cell_Sample) %>%
  summarize(count=n()) %>%
  mutate( 
    n=1,
    mean=count,
    sd=NA,
    se=NA
  ) %>%
  select(-count) 

tx<-  dt %>%
  filter(grepl("Treg",cell_subset ) ) %>%
  group_by(cell_subset, chain,Mice_ID, Cell_Sample) %>%
  summarize(count=n()) %>%
  group_by(cell_subset, chain,Cell_Sample) %>%
  summarise( 
    n=n(),
    mean=mean(count),
    sd=sd(count)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  select(cell_subset,chain,Cell_Sample,n,mean,sd,se) %>%
  mutate(Cell_Sample=ifelse(Cell_Sample=="indiv", "Mouse","Metamouse"))

##a
all_t<-rbind(t1,t2,t3,tx)
all_t$cell_subset<- gsub('s', '', all_t$cell_subset)

all_t %>%
  mutate(cell_subset=factor(cell_subset, levels=c("nTreg","amTreg"))) %>%
  mutate(Cell_Sample=factor(Cell_Sample, levels=c( "Mouse","Metamouse","MetaTCR 1","MetaTCR 2","MetaTCR 3"))) %>%
  mutate(group=ifelse(grepl("MetaTCR", Cell_Sample), "MetaTCR",Cell_Sample)) %>%
  
  ggplot(aes(x=Cell_Sample, y=mean, fill=group, pattern=Cell_Sample)) +
  geom_bar_pattern( stat="identity",
                    color = "black", 
                    alpha=.9,
                    pattern_fill = "black",
                    pattern_angle = 45,
                    pattern_density = 0.1,
                    pattern_spacing = 0.025,
                    pattern_key_scale_factor = 0.6) + 
  scale_pattern_manual(values =  c("none", "none", "stripe", "circle", "crosshatch")) +
  geom_errorbar( aes(x=Cell_Sample, ymin=mean-sd, ymax=mean+sd),
                 width=0.2, colour="black", alpha=0.9, linewidth=.5)+
  facet_grid(chain~cell_subset)+
  scale_fill_manual(values=nejm_color)+
  theme_Publication()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "right",
        legend.direction = "vertical")+
  xlab("")+ylab("Nb of clonotypes")


##b
dt_meta<- all_metaTCR_aggr_nt
dt_all<- dt_all_nt


dt_all<- dt_all %>% 
  dplyr::select(-targetSequences,-TRV,-TRJ,-cdr3aa,-length_cdr3aa,-VJ,-nb_clonotypes,-nb_sequences,
                -Sample_ID_Std,-cell_number,-sequences_per_clonotypes,-Clin_Group,-sex) %>%
  mutate(Cell_Sample=ifelse(Cell_Sample=="indiv", "Mouse","Metamouse"))

dt<- rbind(dt_meta ,dt_all)
temp <- dt %>% mutate(group2=paste(filename, Cell_Sample,cell_subset, chain, sep="%"))
temp<- temp %>%
  mutate(clonotypes=ifelse(filename=="MetaTCR 3", paste0(clonotypes, Mice_ID), clonotypes))
temp <- data_modelling(temp, row = "group2", col = "clonotypes", value = "cloneCount")

mod <- vegan::renyi(temp, scales = c(1))

mod2<- mod %>% 
  data.frame() %>%
  tibble::rownames_to_column() %>%
  reshape2::melt() %>%
  dplyr::select(-variable) %>%
  tidyr::separate(., rowname, c("filename", "Cell_Sample", "cell_subset", "chain"),sep = "%") %>%
  group_by(cell_subset, chain, Cell_Sample) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate(group=ifelse(grepl("MetaTCR", Cell_Sample), "MetaTCR",Cell_Sample)) %>%
  mutate(cell_subset=factor(cell_subset, levels=c("nTreg","amTreg"))) %>%
  mutate(group=factor(group, levels=c("Mouse","Metamouse","MetaTCR"))) %>%
  mutate(Cell_Sample=factor(Cell_Sample, levels=c( "Mouse","Metamouse","MetaTCR 1","MetaTCR 2","MetaTCR 3"))) 
  
ggplot(mod2, aes(x=Cell_Sample, y=mean, fill=group, pattern=Cell_Sample)) +
geom_bar_pattern( stat="identity",
                  color = "black", 
                  alpha=.9,
                  pattern_fill = "black",
                  pattern_angle = 45,
                  pattern_density = 0.1,
                  pattern_spacing = 0.025,
                  pattern_key_scale_factor = 0.6) + 
scale_pattern_manual(values =  c("none", "none", "stripe", "circle", "crosshatch")) +
geom_errorbar( aes(x=Cell_Sample, ymin=mean-sd, ymax=mean+sd),
               width=0.2, colour="black", alpha=0.9, linewidth=.5)+
facet_grid(chain~cell_subset)+
scale_fill_manual(values=nejm_color)+
theme_Publication()+
theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "right",
      legend.direction = "vertical")+
xlab("")+ylab("Shannon")+
coord_cartesian(ylim=c(7,12))


#Figure S4####
##a
# dt_meta<- readRDS("./files/all_metaTCR_aggr_nt.rds")
# dt_all<- readRDS("./files/dt_all_nt.rds")
dt_all<- dt_all %>% 
  dplyr::select(-targetSequences,-TRV,-TRJ,-cdr3aa,-length_cdr3aa,-VJ,-nb_clonotypes,-nb_sequences,
                -Sample_ID_Std,-cell_number,-sequences_per_clonotypes,-Clin_Group,-sex) %>%
  mutate(Cell_Sample=ifelse(Cell_Sample=="indiv", "Mouse","Metamouse"))

X_CL<- rbind(dt_meta ,dt_all)

X_CL <- X_CL %>%
  group_by(filename, chain,cell_subset) %>%
  arrange(desc(cloneCount)) %>%
  mutate(clone_ID = row_number()) %>%
  mutate(quartile = cut(clone_ID, breaks = 4, labels = paste0("Q", 1:4)))

summary_quartile <- X_CL %>%
  group_by(filename, chain, Cell_Sample, cell_subset, quartile) %>%
  summarise(n_distinct(clonotypes))

X_CL2 <- X_CL %>%
  mutate(group2=paste(filename, Cell_Sample,cell_subset, chain, sep="%")) %>%
  named_group_split(chain , cell_subset,sep="_") 

X_CL2 <- lapply(X_CL2, dkast, gene = "clonotypes", value = "cloneCount")

MH_res2 <- mclapply(X_CL2, Compute_MH, method = "ji", out="df", mc.cores = 4) ## method mh or ji

MH_res3<-rbindlist(MH_res2, idcol = "group")

MH_res3 <- MH_res3 %>%
  separate(., Var1, into = c( "filename1", "Cell_Sample1", "cell_subset1","chain1"), sep = "%") %>%
  separate(., Var2, into = c( "filename2","Cell_Sample2", "cell_subset2","chain2"), sep = "%") %>%
  separate(., group, into = c("chain", "cell_subset"), sep = "_")

rm<- c("MetaTCR10","MetaTCR2","MetaTCR3","MetaTCR4","MetaTCR5","MetaTCR6","MetaTCR7","MetaTCR8", "MetaTCR9" )

MH_res4 = MH_res3 %>%
  filter(!grepl(paste0(rm, collapse="|"), filename1)) %>%
  filter(!grepl(paste0(rm, collapse="|"), filename2)) %>%
  filter(Cell_Sample1 != Cell_Sample2) %>% #### If we don't want to compare indiv sample with pool samples
  mutate(keep=ifelse(Cell_Sample1=="Mouse" | Cell_Sample2=="Mouse", "keep","no" )) %>%
  filter(keep=="keep") %>%
  mutate(group=ifelse(Cell_Sample1=="Mouse" , Cell_Sample2, Cell_Sample1)) %>%
  mutate(cell_subset=factor(cell_subset, levels=c("nTreg","amTreg"))) %>%
  mutate(group2=ifelse(grepl("MetaTCR", group), "MetaTCR",group)) %>%
  select(chain, cell_subset, value, group, group2) %>%
  distinct()


stat.test <- MH_res4 %>%
  group_by(chain, cell_subset) %>%
  rstatix::wilcox_test(value ~ group) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_y_position()

MH_res4 %>%
  ggplot(.,  aes(x=group, y=value)) +
  geom_boxplot_pattern( aes(fill=group2, pattern=group),
                        outlier.shape = NA,
                        linewidth=.3,
                        color = "black", 
                        alpha=.9,
                        pattern_fill = "black",
                        pattern_angle = 45,
                        pattern_density = 0.1,
                        pattern_spacing = 0.025,
                        pattern_key_scale_factor = 0.6) + 
  scale_pattern_manual(values =  c("none", "stripe", "circle", "crosshatch"))+
  geom_jitter(fill="black",shape=21,size=1,width=.2,alpha=.5, color="black")+
  facet_grid(chain~cell_subset) +
  scale_fill_manual(values= c("#20854e","#bc3c29"))+
  stat_pvalue_manual(data = stat.test, label = "p.adj.signif", tip.length = 0.01 ,  hide.ns = TRUE) +
  theme_Publication()+
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5))+
  xlab("")+
  ylab("Jaccard")

 
##b
 
 X_CL_Q1 <- X_CL %>% filter(quartile == "Q1") #%>%  filter(filename!="tripod-64-1932_R1.txt") 
 X_CL_Q1 <- X_CL_Q1 %>% 
   mutate(group2=paste(filename, Cell_Sample,cell_subset, chain, sep="%")) %>%
   named_group_split(chain,cell_subset,  sep="%") 
 X_CL_Q1 <- lapply(X_CL_Q1, dkast, gene = "clonotypes", value = "freq_clo")
 
 data_cor <- lapply(X_CL_Q1, function(x) x <- cor(x))

 corr_all <- lapply(data_cor,   melt)
 corr_all<-rbindlist(corr_all, idcol = "group")
 
 corr_all <- corr_all %>%
   separate(., Var1, into = c( "filename1", "Cell_Sample1", "cell_subset1","chain1"), sep = "%") %>%
   separate(., Var2, into = c( "filename2","Cell_Sample2", "cell_subset2","chain2"), sep = "%") %>%
   separate(., group, into = c("chain", "cell_subset"), sep = "%")
 
 rm<- c("MetaTCR10","MetaTCR2","MetaTCR3","MetaTCR4","MetaTCR5","MetaTCR6","MetaTCR7","MetaTCR8", "MetaTCR9" )
 corr_all2 = corr_all %>%
   filter(!grepl(paste0(rm, collapse="|"), filename1)) %>%
   filter(!grepl(paste0(rm, collapse="|"), filename2)) %>%
   filter(Cell_Sample1 != Cell_Sample2) %>% #### If we don't want to compare indiv sample with pool samples
   mutate(keep=ifelse(Cell_Sample1=="Mouse" | Cell_Sample2=="Mouse", "keep","no" )) %>%
   filter(keep=="keep") %>%
   mutate(group=ifelse(Cell_Sample1=="Mouse" , Cell_Sample2, Cell_Sample1)) %>%
   mutate(cell_subset=factor(cell_subset, levels=c("nTreg","amTreg"))) %>%
   mutate(group2=ifelse(grepl("MetaTCR", group), "MetaTCR",group)) %>%
   select(chain, cell_subset, value, group, group2) %>%
   distinct()
 

 stat.test <-corr_all2 %>%
   group_by(chain, cell_subset) %>%
   rstatix::wilcox_test(value ~ group) %>%
   rstatix::adjust_pvalue(method = "holm") %>%
   rstatix::add_significance("p.adj") %>%
   rstatix::add_y_position()
 
 
corr_all2 %>%
   ggplot(.,  aes(x=group, y=value)) +
   geom_boxplot_pattern( aes(fill=group2, pattern=group),
                         outlier.shape = NA,
                         color = "black", 
                         linewidth=.3,
                         alpha=.9,
                         pattern_fill = "black",
                         pattern_angle = 45,
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6) + 
   scale_pattern_manual(values =  c("none", "stripe", "circle", "crosshatch"))+
   geom_jitter(fill="black",shape=21,size=1,width=.2,alpha=.5, color="black")+
   facet_grid(chain~cell_subset) +
   scale_fill_manual(values= c("#20854e","#bc3c29"))+
   stat_pvalue_manual(data = stat.test, label = "p.adj.signif", tip.length = 0.01 , hide.ns = TRUE) +
   theme_Publication()+
   theme(legend.position = "right", legend.direction = "vertical",
         axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5))+
   xlab("")

#Figure 3####

##a
for(var_chain in c("TRA", "TRB")) {
  for (i in unique(dt_filtered$filename)) {
    print(i)
    newdata <- dt_filtered %>% 
      filter(filename == i & chain == var_chain) %>%
      dplyr::group_by(Mice_ID, cell_subset, Cell_Sample, chain, cdr3aa, filename) %>%
      summarise(freq_clo = sum(freq_clo)) %>%
      ungroup()
    
    info <- newdata %>% select(Mice_ID, cell_subset, Cell_Sample, chain) %>% distinct()
    
    newdata <- newdata %>% 
      dplyr::slice_max(n = 1000, order_by = freq_clo, with_ties = F) %>%
      dplyr::select(cdr3aa, freq_clo)
    
    g1s <- find_pairs(newdata$cdr3aa, newdata$cdr3aa)
    
    g1s <- g1s %>%
      graph_from_data_frame(directed = F) %>% simplify()
    
    ## Create a plot  
    V(g1s)$degree <- degree(g1s)
    cc <- clusters(g1s)
    V(g1s)$community <- cc$membership
    V(g1s)$freq <-newdata$freq_clo
    
    # 2. Create a plot
    ggraph(g1s, layout = "graphopt") +
            geom_edge_link(alpha = 0.7) + 
            geom_node_point(aes(fill = degree, size =freq ), alpha = 0.8, shape = 21, show.legend = T, color = "black") +
            #scale_color_distiller(guide = F, palette = "Paired") +
            scale_fill_gradient(low="blue", high="red", limits=c(0,40))+
            scale_size_continuous(limits =c(0, 1))+
            ggtitle(paste(info$Mice_ID, info$cell_subset, info$Cell_Sample, info$chain, sep = " "))+
            theme_graph()
    cc <- clusters(g1s)
    cc$csize <- cc$csize[cc$membership]
    networks.stats_2b <- data.frame(cdr3aa = names(V(g1s))) %>%
      merge(data.frame(cdr3aa = names(cc$membership),
                       cluster_id = cc$membership,
                       degree = degree(g1s, mode = "all"),
                       cluster_size = cc$csize,
                       Clustering_coef_each = transitivity(g1s, type = "localundirected", isolates = "zero"),
                       Clustering_coef_full = transitivity(g1s, type = "global", isolates = "zero")))

    newdata <- merge(newdata, networks.stats_2b, by = "cdr3aa", all.x = T, all.y = F) %>% distinct(cdr3aa, freq_clo, cluster_id, degree, cluster_size, Clustering_coef_full, Clustering_coef_each)
    write.table(newdata, paste0("./files/",i, "_CDR3aa_info_", var_chain, ".txt"), sep = "\t", row.names = F, quote = F)
  }
}

## b
# file_list <- list.files(path = "./files/", pattern = "_CDR3aa_info_", include.dirs = T)

networks_filename = file_list
names(file_list) <- file_list
# file_list <- lapply(paste0("./files/", file_list), fread)

file_list <- mapply(cbind, file_list, "filename" = networks_filename, SIMPLIFY = F)
networks_dt = rbindlist(file_list)

networks_dt = separate(networks_dt, filename, into = c("filename", "chain"), sep = "_CDR3aa_info_")
networks_dt = separate(networks_dt, chain, into = c("chain"), sep = "[.]")

networks_dt = merge(networks_dt, unique(metadata[,-c("chain", "cell_number", "Clin_Group", "sex", "group2")]), by = "filename")

## Graph 2
cluster_plot <- networks_dt %>%
  dplyr::group_by(filename, chain) %>% 
  arrange(filename, Clustering_coef_each, chain) %>%
  mutate(rank = rank(-Clustering_coef_each, ties.method =  "first")) %>%
  mutate(AUC = MESS::auc(rank, Clustering_coef_each, type = 'spline')) #%>%

cluster_plot$Cell_Sample <- factor(cluster_plot$Cell_Sample, levels = c("indiv", 'pool'))
cluster_plot$cell_subset <- factor(cluster_plot$cell_subset, levels = c("CD8+", 'CD4+', 'nTreg', "amTreg"))

stat.test <- cluster_plot %>% select(Cell_Sample, AUC, chain, cell_subset) %>% distinct() %>%
  group_by(cell_subset, chain) %>%
  rstatix::wilcox_test(AUC ~ Cell_Sample) %>%
  #ggpubr::compare_means(data = MH_res2, group.by = c("cell_subset1", "chain1"), formula =  value ~ Cell_Sample1, method = "wilcox") %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>% rstatix::add_xy_position()


cluster_plot %>% select(Cell_Sample, AUC, chain, cell_subset) %>% distinct() %>%
  ggpubr::ggboxplot(x = "Cell_Sample", y = "AUC", alpha=.5,  fill = "Cell_Sample", color = "Cell_Sample", outlier.shape = NA,
                    facet.by = c("chain","cell_subset"), scales ="free_y")+
  geom_jitter( aes( fill=Cell_Sample),size=1.5,color = "black", shape=21,alpha=.4, position = position_jitterdodge(.2))+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  ylab("Clustering coeffecient (AUC)") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01)


##c

networks_dt3 <- networks_dt %>% dplyr::group_by(chain, cell_subset, Cell_Sample) %>%
  dplyr::mutate(total_mice = n_distinct(Mice_ID)) %>%
  dplyr::group_by(cdr3aa, chain, cell_subset, Cell_Sample) %>%
  dplyr::mutate(sharing = n_distinct(Mice_ID)) %>%
  dplyr::mutate(cdr3aa_category = if_else(sharing == 1 , "private", "public")) %>%
  dplyr::group_by(chain, cell_subset, Cell_Sample, cdr3aa_category, Mice_ID) %>%
  dplyr::summarise(mean_Clustering_coef_each = mean(Clustering_coef_each))


## Recodage de networks_dt3$Cell_Sample
networks_dt3$Cell_Sample[networks_dt3$Cell_Sample == "indiv"] <- "Mouse"
networks_dt3$Cell_Sample[networks_dt3$Cell_Sample == "pool"] <- "Metamouse"

## Réordonnancement de networks_dt3$cell_subset
networks_dt3$cell_subset <- factor(networks_dt3$cell_subset,
                                   levels = c("CD8+", "CD4+", "nTreg", "amTreg")
)

## Réordonnancement de networks_dt3$cell_subset
networks_dt3$Cell_Sample <- factor(networks_dt3$Cell_Sample,
                                   levels = c("Mouse", "Metamouse")
)


stat.test <- networks_dt3 %>% 
  group_by(cell_subset, chain, cdr3aa_category) %>%
  rstatix::wilcox_test(mean_Clustering_coef_each ~ Cell_Sample) %>%
  #ggpubr::compare_means(data = MH_res2, group.by = c("cell_subset1", "chain1"), formula =  value ~ Cell_Sample1, method = "wilcox") %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>% rstatix::add_xy_position(x = "cdr3aa_category", dodge = 0.8)


networks_dt3 %>%
  ggpubr::ggboxplot(x = "cdr3aa_category", y = "mean_Clustering_coef_each", alpha=.5,  fill = "Cell_Sample", color = "Cell_Sample", outlier.shape = NA,
                    facet.by = c("chain","cell_subset"), scales ="free_y")+
  geom_jitter(aes(fill=Cell_Sample),size=1.5,color = "black", shape=21,alpha=.4, position = position_jitterdodge(.1))+
  scale_fill_manual(values = nejm_color)+
  scale_color_manual(values = nejm_color)+
  theme_Publication()+
  ylab("Mean clustering coefficient") +
  facet_grid(chain~cell_subset, scales = "free")+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggpubr::stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01)








