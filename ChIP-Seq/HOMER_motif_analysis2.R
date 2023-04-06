#########################
# Visualization of HOMER output
# Fig 3 c,d, S3 b



### density plot
library(tidyverse)

meis2_all_density - read.table(.HOMERmotif_density_meis2peaks.txt, skip=1,header=F, sep=t)

columns - c(distance , hexa_total_sites	,hexa_pos_sites, hexa_neg_sites,	deca_total_sites,	deca_pos_sites, deca_sites,
             	A_frequency, C_frequency, G_frequency,	T_frequency)

names(meis2_all_density) - columns


meis2_motif_enrichment2 - ggplot(meis2_all_density, aes(x=distance)) +
  geom_smooth(aes(y=hexa_total_sites), se=F, method=gam, color=#F47E70)+
  geom_smooth(aes(y=deca_total_sites), se=F, method=gam, color=#8FD1C4) +
  theme_classic()+
  labs(x=distance to peak summit [bp], y=relative motif enrichment) +
  theme(axis.text = element_text(size=7), axis.title=element_text(size=8))
meis2_motif_enrichment2




### motif enrichment plot
library(tidyverse)
library(RColorBrewer)

#### results from HOMER


############ all MEIS2 peaks

homer_enrichment_meis2_all - read.table(.HOMERknownResults_MEIS2_peaks_all.txt, 
                                   header=F,row.names=NULL,sep=t, skip=1) 

homer_enrichment_meis2_all[1,4]
names(homer_enrichment_meis2_all) - c(motif_name, consensus, p_value, log_pvalue, q_value_Benjamini, no_of_targets_seqs_with_motif, percent_of_target_seqs_with_motif, no_of_background_seqs_with_motif, percent_of_background_seqs_with_motif)

#convert ln value to log10
homer_enrichment_meis2_all - homer_enrichment_meis2_all %%
  mutate(n_log10_pval = -round((log_pvalue  log10(2.71828) ),1))
  
# "log p-value" column is natural logarithm!

homer_enrichment_all2 - homer_enrichment_meis2_all %%
  mutate(percent_of_target_seqs_with_motif = as.numeric(gsub(%,,percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub(%,,percent_of_background_seqs_with_motif))) %%
  mutate(fold_enrichment = percent_of_target_seqs_with_motifpercent_of_background_seqs_with_motif) %%
  tidyrseparate(motif_name,c( motif_short,motif_info1, motif_info2), sep=(, remove=F) %%
  dplyrselect(-c(motif_info1,motif_info2)) %%
  mutate(rownumber = seq(1nrow(.)))


# manually selected sequences from HOMER knownResults, to avoid redundancy
rowselectionB - c(1,4, 7,8,12, 16,20,13,69,15 )
homer_enrichment_all_selectedB - homer_enrichment_all2 %%
  filter(rownumber %in% rowselectionB) %%
  mutate() %%
  arrange(factor(rownumber, levels = rev(rowselectionB))) %%
  mutate(motif_short = as.factor(motif_short)) %%
  mutate(motif_short = factor(motif_short, 
                              levels=c(Nanog , Sp5 ,   Tbx5  , Oct6 ,  Isl1  , Nkx6.1 ,Lhx2  , Dlx3 ,  Pbx3  , Meis1 )))

# homer_enrichment_plot 
homer_enrichment_plot_selectionB - ggplot(homer_enrichment_all_selectedB, aes(x= motif_short))+
  theme_bw()+
  geom_col(aes(y=percent_of_background_seqs_with_motif),fill=orange,alpha=0.5, color=red)+
  geom_col(aes(y=percent_of_target_seqs_with_motif),fill=lightgrey, alpha=0.5, color=black)+
  theme(axis.line = element_line(colour = black),panel.grid = element_blank(), panel.border =  element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(face=bold, size=9), axis.text.x= element_text(size=9), axis.title.x=element_text(size=10))+
  geom_text(aes(y=percent_of_target_seqs_with_motif, label=paste((,round(fold_enrichment,2),), sep=),hjust=-0.1),size=90.36)+
  ylim(0,100)+
  labs(x=NULL, y=percent of peaks with motif)+
  coord_flip()
homer_enrichment_plot_selectionB


###########################################
############ divide by enhancers and promoters


############# MEIS2 peaks in promoters 


homer_enrichment_promoter - read.table(.HOMERknownResults_MEIS2_peaks_promoter.txt, 
                                   header=F,row.names=NULL,sep=t, skip=1)
names(homer_enrichment_promoter) - c(motif_name, consensus, p_value, log_pvalue, q_value_Benjamini, no_of_targets_seqs_with_motif, percent_of_target_seqs_with_motif, no_of_background_seqs_with_motif, percent_of_background_seqs_with_motif)

homer_enrichment_promoter2 - homer_enrichment_promoter %%
  mutate(percent_of_target_seqs_with_motif = as.numeric(gsub(%,,percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub(%,,percent_of_background_seqs_with_motif))) %%
  mutate(fold_enrichment = percent_of_target_seqs_with_motifpercent_of_background_seqs_with_motif) %%
  tidyrseparate(motif_name,c( motif_short,motif_info1, motif_info2), sep=(, remove=F) %%
  dplyrselect(-c(motif_info1,motif_info2)) %%
  mutate(rownumber = seq(1nrow(.)))

# log p-value is natural logarithm!

rowselection_promoter - c(1,3,7,8,9,13,15,17,18,27)
homer_enrichment_promoter_selected - homer_enrichment_promoter2 %%
  filter(rownumber %in% rowselection_promoter) %%
  mutate()



##############  MEIS2 peaks in ennhancers 

homer_enrichment_enhancer - read.table(homevolker.kittkeDLX_and_MEIS_binding_mouse_embryooutputmotif_analysisHOMERMEIS2_200bp_enhancer_peaks_newknownResults.txt, 
                                        header=F,row.names=NULL,sep=t, skip=1)
names(homer_enrichment_enhancer) - c(motif_name, consensus, p_value, log_pvalue, q_value_Benjamini, no_of_targets_seqs_with_motif, percent_of_target_seqs_with_motif, no_of_background_seqs_with_motif, percent_of_background_seqs_with_motif)

homer_enrichment_enhancer2 - homer_enrichment_enhancer %%
  mutate(percent_of_target_seqs_with_motif = as.numeric(gsub(%,,percent_of_target_seqs_with_motif)), percent_of_background_seqs_with_motif = as.numeric(gsub(%,,percent_of_background_seqs_with_motif))) %%
  mutate(fold_enrichment = percent_of_target_seqs_with_motifpercent_of_background_seqs_with_motif) %%
  tidyrseparate(motif_name,c( motif_short,motif_info1, motif_info2), sep=(, remove=F) %%
  dplyrselect(-c(motif_info1,motif_info2)) %%
  mutate(rownumber = seq(1nrow(.)))


# enhancers-manually selected sequences
rowselection_enhancer - c(2,4,6,7,10,13,15,16,17,19)
homer_enrichment_enhancer_selected - homer_enrichment_enhancer2 %%
  filter(rownumber %in% rowselection_enhancer) %%
  mutate()


##################################
##### grouped barplot for easier comparison between enhancers and promoters
# 1) rename columns
# 2) bind rows to generate one table, add ID
# 3) 



homer_enrichment_promoter_enhancer_combined - homer_enrichment_promoter2 %%
  bind_rows(homer_enrichment_enhancer2,.id = ID) %%
  mutate(ID = ifelse(ID == 1, promoter, enhancer))


# select motifs by name

motifs - c(Meis1, Pbx3, Dlx3, Lhx2, Nkx6.1, Isl1, Oct6, Tbx5, Sp5, Nanog)


homer_enrichment_promoter_enhancer_selected - homer_enrichment_promoter_enhancer_combined %%
  filter(motif_short %in% motifs)

# convert to wide
homer_enrichment_promoter_enhancer_wide - homer_enrichment_promoter_enhancer_selected %%
  dplyrselect(c(14,12)) %%
  pivot_wider(names_from = ID, values_from = fold_enrichment) %%
  mutate(promoterVSenhancer = enhancerpromoter)



homer_enrichment_plot_promoter_enhancer_selected - ggplot(homer_enrichment_promoter_enhancer_selected, 
                                                           aes(x= factor(motif_short, levels=rev(motifs)), fill= ID))+
  theme_bw()+
  geom_col(aes(y=percent_of_background_seqs_with_motif ),width=0.8,alpha=1, position=dodge, color=red)+
  geom_col(aes(y=percent_of_target_seqs_with_motif),width = 0.8, alpha=0.5, position=dodge, color=black)+
  geom_text(aes(y=percent_of_target_seqs_with_motif, label=paste((,round(fold_enrichment,2),), sep=)),size=80.36, hjust=-0.1, 
            position=position_dodge(width=0.8))+
  ylim(0,100)+
  labs(x=NULL, y=% of peaks with motif, fill=peak category)+
  scale_fill_manual(values = c(#FDBF6F, #CAB2D6))+
  theme(axis.line = element_line(colour = black),panel.grid = element_blank(), panel.border =  element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(face=bold, size=9), axis.text.x= element_text(size=9), axis.title.x=element_text(size=9),
        legend.position = right, legend.direction = vertical)+
  coord_flip()
homer_enrichment_plot_promoter_enhancer_selected