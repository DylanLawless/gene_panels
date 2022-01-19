################################################
####  Install or load packages              ####
################################################

library(dplyr)
library(stringr)
library(ggplot2)
require(gridExtra)
library(grid)
library(ggsignif)


################################################
####  Import the dataset and attach         ####
################################################

# Dataset: 
gnomad_pLI <- read.csv(file = "./gnomAD_pLI.sort.csv",sep=",",header=T)
immune <- read.csv(file = "./immunology_genes_symbols.csv",sep=",",header=T)

df <- left_join(immune, gnomad_pLI,
          by = "gene", 
          copy=FALSE
          ) 



left_join(x, y, by = NULL, copy=FALSE, suffix=c(“.x”,“.y”),...) Join matching values from y to x.









mrf.present <- mrf[which(mrf[,"Boolean_conservation_score"]<1),]

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

################################################
####  Part (i) MRF score comparison         ####
################################################

p1 <-
  completeFun(mrf.present, "Variant_reported") %>%
  ggplot(aes(x=reported_pathogenic, y=Mutation_rate_residue_frequency))+
  scale_x_discrete(labels=c("reported" = "Known\npathogenic ", 
                            "unreported" = "Unreported"))+
  geom_boxplot(alpha=0)+
  geom_point(aes(fill = reported_pathogenic:Gene),
             size=2, shape=21, 
             position = position_jitter(width = 0.2, height = 0.0005))+
  scale_fill_manual(values=c("#FC600A", "#195FE0", "#8f7e76", "#707989"))+ 
  facet_grid(~Gene)+
  labs(x = "Status as a pathogenic variant",
       y = "MRF score")+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect("#F7F7F7"))+
  geom_signif(comparisons = list(c("reported", "unreported")), 
              map_signif_level=TRUE)+
  ggtitle("(i)")

################################################
####  Part (ii) allele frequency comparison ####
################################################

p2 <-
  popgen_patho %>%
  ggplot(aes(x=reported_pathogenic, y=AF))+
  scale_x_discrete(labels=c("reported" = "Known\npathogenic ", 
                            "unreported" = "Unreported"))+
  geom_boxplot(alpha=0)+
  geom_point(aes(fill = reported_pathogenic:Gene),
             size=2, shape=21, 
             position = position_jitter(width = 0.2, height = 0.0000004))+
  scale_fill_manual(values=c("#FC600A", "#195FE0", "#8f7e76", "#707989"))+
  facet_grid(~Gene)+
  labs(x = "Status as a pathogenic variant",
       y = "Allele frequency")+
  theme(legend.position="none",
        axis.title.x=element_blank(), 
        panel.background = element_rect("#F7F7F7"))+
  geom_signif(comparisons = list(c("reported", "unreported")), 
              map_signif_level=TRUE)+
  ggtitle("(ii)")

################################################
####      Arrange the plots together        ####
################################################

grid.arrange(p1, p2, ncol=2, bottom = textGrob("Status as a pathogenic variant and present in general population"))
# size 5 x 8

################################################
####  Figure 3                              ####
####  "RAG1 and RAG2 MRF score categories   ####
####  and variants assayed to date.         ####
################################################
####  To show how many of the top mrf score ####
####  candidates have been tested in vivo   ####
####  or in vitro, and show by colour       ####
####  their recombination activity (%)      ####
################################################

# Dataset: The percentage of recombination activity (SEM) per residue mutation tested in vitro by Tirosh et al 2018, Lee et al 2014, Lawless et al 2018.
activity <- read.csv(file = "./recombination_activity.csv",sep=",",header=T) # for activity, where a single residue was tested with alternative mutations the median was used.

# Count mrf categories
mrf.count <-  mrf %>% 
  dplyr::group_by(Mutation_rate_residue_frequency, Gene) %>%
  dplyr::summarise(length(Mutation_rate_residue_frequency))

# Add the functionally validated recombination activiy data
mrf.activity <-  left_join(mrf, activity)

# Count the number of activity scores per mrf category
mrf.activity.summary <-
  mrf.activity %>% 
  dplyr::group_by(Mutation_rate_residue_frequency, Gene, activity) %>%
  dplyr::summarise(count =length(activity))

# Round the decimal places to 4 digits for graphing on x-axis
mrf.activity.summary[, "score_rounded"] <- 
  round(mrf.activity.summary[,"Mutation_rate_residue_frequency"], digits=4)

# Graph the activity measured per mrf category group. Untested variants are shown in grey.
ggplot() + 
  geom_bar(data = mrf.activity.summary %>%  
             # convert the MRF score labels fron values to characters to remove x-axis gaps
             dplyr::mutate(score_char = as.character(score_rounded)), 
           aes(score_char, count, fill = activity, colour = ""), 
           width = 1, stat = "identity", position = "stack") + 
  # plot the same bars with no fill (alpha=0), but colour black to provide bar boarders
  geom_bar(data = mrf.activity.summary %>%  dplyr::mutate(score_char = as.character(score_rounded)), aes(score_char, count), width = 1, stat = "identity", position = "stack",
           colour="black", alpha=0) + 
  # set colour scale with grey value for untested residues (NA) and rename legends
  scale_fill_gradient(low="brown3", high="cornsilk1", 
                      na.value = "grey",
                      name = "Recombination\nactivity (%)") + 
  scale_colour_manual(values=NA) +              
  guides(colour=guide_legend("Untested", override.aes=list(colour="grey")))+
  # theme and labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect("#F7F7F7"))+
  facet_grid(~Gene, scales = "free")+
  labs(x ="MRF score categories \n(Likelihood of presentation)",
       y = "Number of residues")+
  annotate("text", x = 8, y = 120, label = expression(decreasing %<->% increasing))

################################################
####  Figure 4                              ####
####  "False positives in Transib domains   ####
####  do not worsen probability prediction" ####
################################################
####  To present MRF score for residues in  ####
####  the Transib region and seperate into  ####
####  their three categories; show the % of ####
####  each category mutated on GnomAD;      ####
####  lastly apply conservation score C to  ####
####  illustrate that simulated false       ####
####  positives in even the most critical   ####
####  motifs are non-significant when the   ####
####  individual conservation rates per     ####
####  residue are accounted for             ####
################################################

# Dataset: Transib core motifs with mrf scores per category.
transib <- read.csv(file = "./transib.csv",sep=",",header=T)

# Dataset: Percentage mutated per group in general population (GnomAD).
transib_pc <- read.csv(file = "./transib_pc.csv",sep=",",header=T)

trasibItalic <- expression(paste(italic("Transib"), "motifs")) # this sets the italics for the figure label.

p3 <- 
  transib %>%
  ggplot(aes(x=status, y=mrf))+
  geom_boxplot(alpha=0.2)+
  geom_point(aes(fill = status), 
             size=3, shape=21, 
             position = position_jitter(width = 0.3, height = 0.0001))+
  scale_fill_manual(values=c("#006080", "#39acac", "#d9e4cd"))+  
  labs(x=trasibItalic, y = "MRF score")+
  theme(legend.position="none", panel.background = element_rect("#F7F7F7"))+
  geom_signif(stat="identity",
              data=data.frame(x=c(3, 1), xend=c(1, 2),
                              y=c(.045, .048), annotation=c("***", "*")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  geom_signif(comparisons=list(c("S1", "S2")), annotations="***",
              y_position = 9.3, tip_length = 0, vjust=0.4)+
  ggtitle("(i)")+
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))

p4 <- 
  transib_pc %>%
  ggplot(aes(x=status, y=mrf))+
  geom_bar(aes(fill=status), stat="identity", alpha=0.8, colour = "black")+
  scale_fill_manual(values=c("#006080", "#39acac", "#d9e4cd"))+  
  labs(x=trasibItalic, y = "Residues mutated \non GnomAD (%)")+
  theme(legend.position="none", panel.background = element_rect("#F7F7F7"))+
  ggtitle("(ii)")+
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))

p5 <- 
  transib[which(transib[,"mrf"]>0),] %>%
  ggplot(aes(x=status, y=mrf))+
  geom_boxplot(alpha=0.2)+
  geom_point(aes(fill = status), 
             size=3, shape=21, 
             position = position_jitter(width = 0.3, height = 0.0001))+
  scale_fill_manual(values=c("#006080", "#39acac", "#d9e4cd"))+  
  labs(x=trasibItalic, y = "MRF score")+
  theme(legend.position="none", panel.background = element_rect("#F7F7F7"))+
  geom_signif(stat="identity",
              data=data.frame(x=c(3), xend=c(1),
                              y=c(.048), annotation=c("ns")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  geom_signif(comparisons=list(c("S1", "S2")), annotations="***",
              y_position = 9.3, tip_length = 0, vjust=0.4)+
  ggtitle("(iii)")+
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))

################################################
####      Arrange the plots together        ####
################################################
grid.arrange(p3, p4, p5, ncol=3)
# size 6.2 x 12

################################################
####  Figure 5                              ####
################################################
####  "A linear regression model of RAG1/2  ####
####  MRF scoring in cases of primary       ####
####  immune deficiency."                   ####
################################################

# Import the dataset and attach.
RAG_mrf_human_count <- read.csv(file = "./RAG_mrf_human_count.csv",sep=",",header=T)

# Display dataset as a scatterplot with log transformed x axis.
RAG_mrf_human_count %>%
  ggplot(aes(x=log(count), y=mrf, color=mutant, size=count))+
  geom_point(alpha=1)+
  geom_smooth(method = lm, alpha = 0.1, show.legend = FALSE)+
  facet_grid(mutant ~ Protein)+ 
  labs(x = "log(Number of human cases)",
       y = "MRF",
       size = "Number of \nhuman cases",
       color = "Variant in \nPID patient")+
  scale_color_manual(values=c("#FC600A", "#FFD438"))+
  theme(panel.background = element_rect("#F7F7F7"))

#### ----------- Damaging ----------- ####
#### MRF scores in damaging variants ####
# Make dataframe for RAG1 damaging variants only and display the slope.
R1_damaging <- RAG_mrf_human_count %>%
  select(count, mrf, Protein, mutant) %>%
  filter(Protein == "RAG1") %>%
  filter(mutant == "Damaging")
lm(mrf ~ count, data = R1_damaging)

# Make dataframe for RAG2 damaging variants only and display the slope.
R2_damaging <- RAG_mrf_human_count %>%
  select(count, mrf, Protein, mutant) %>%
  filter(Protein == "RAG2") %>%
  filter(mutant == "Damaging")
lm(mrf ~ count, data = R2_damaging)

# Get the summary information for these linear models.
summary(lm(mrf ~ count, data = R1_damaging))
summary(lm(mrf ~ count, data = R2_damaging))

#### ----------- Non-damaging ----------- ####
#### MRF scores in non-damaging variants ####
# Make dataframe for RAG1 damaging variants only and display the slope.
R1_non_damaging <- RAG_mrf_human_count %>%
  select(count, mrf, Protein, mutant) %>%
  filter(Protein == "RAG1") %>%
  filter(mutant == "Non-damaging")
lm(mrf ~ count, data = R1_non_damaging)

# Make dataframe for RAG2 damaging variants only and display the slope.
R2_non_damaging <- RAG_mrf_human_count %>%
  select(count, mrf, Protein, mutant) %>%
  filter(Protein == "RAG2") %>%
  filter(mutant == "Non-damaging")
lm(mrf ~ count, data = R2_non_damaging)

# Get the summary information for these linear models.
summary(lm(mrf ~ count, data = R1_non_damaging))
summary(lm(mrf ~ count, data = R2_non_damaging))

################################################
####  Figure 6 illustrated using GraphPad   ####
####  Prism using PHRED-CADD scores from    ####
####  Kircher et al., the mrf scores with   ####
####  1% window and 75th percentile cutoff. ####
################################################

################################################
####  Figure 7                              ####
####  "RAG1 PHRED-scaled CADD score versus  ####
####  MRF score against HGMD data.          ####
################################################
####  Import HGMD data                      ####
################################################

# Dataset: Variants reported on HGMD with citations, phenotype, accession numbers.
hgmd <- read.csv(file = "./hgmd_variants.csv",sep=",",header=T)

# "HGMD" data for RAG1 and RAG2 sourced from The Human Gene Mutation Database at the Institute of Medical Genetics in Cardiff
# http://www.hgmd.cf.ac.uk/ac/index.php

# We want to compare (i) the mrf scores of variants in the general population which are not repoted as pathogenic
# to (ii) the mrf score of variants reported on hgmd.

# (i) For under the "reported_pathogenic" tab only use unreported 
# We do this first seperately rather than joining the hmgb data since lots of hmgb variant are from the primamary literature sources. 

mrf.unreported <- mrf %>% filter(str_detect(reported_pathogenic, "unreported"))

# (ii) then we add the mrf scores for variants reported in hgmd.

hgmd.mrf <- merge(mrf, hgmd)
hgmd.mrf.1 <- hgmd.mrf[which(hgmd.mrf[,"Boolean_conservation_score"]<1),]

# The two sets now have their mrf score and we can combine; variants in the general population that are not reported as pathogenic and variants reported on hgmd.

joined <- hgmd.mrf.1 %>% full_join(mrf.unreported)

# hgmd reports multiple alt variants for the same residues, these need only be counted once since the mrf score is set at the residue level and would change outcome if counted twice.

joined.dedup <- joined[!duplicated(joined[c("Gene","Residue","hgmd")]),]

# Shorten the dataset to read more easily.

joined.dedup.short <- joined.dedup[,c("Gene", "Residue", "Mutation_rate_residue_frequency", "hgmd", "Boolean_conservation_score", "reported_pathogenic")]

# After full join, the mrf data for variants not reported in hgmd have "NA" in the hgmd column, these are unreported and should be labelled thusly (reported variants are labelled "1")

joined.dedup.short$hgmd[is.na(joined.dedup.short$hgmd)] <- "unreported"

################################################
####  Figure 7                              ####
####  "RAG1 PHRED-scaled CADD score versus  ####
####  MRF score against HGMD data.          ####
################################################
####  Import CADD data                      ####
################################################

# CADD scores are calculated per nucleotide. Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892.
# Data was sourced from "All possible SNVs" from whole genome data then extracting the data for coding regions of RAG1 and RAG2
# https://cadd.gs.washington.edu/download

# Most people would like a tool that predicts variant pathogenicity.
# Since MRF is focused on variant probability, it is useful to combine MRF scores with CADD score, a commonly used for "predicting the deleteriousness of variants".
# CADD scores are based on nucleotide level. We are interested in approximate cadd scores per codon.
# For every nucleotide position there are 3 alternative variants to consider; e.g.
#Chrom	Pos	Ref	Alt1	Alt2	Alt3	RawScore1	RawScore2	RawScore3	PHRED1	PHRED2	PHRED3
#11	36594855	A	C	G	T	3.021877	2.398474	3.048749	22.3	18.81	22.4
# Both the raw cadd score and PHRED-scaled scores are shown here (as they appear in the raw data from doi: 10.1038/ng.2892)

# To produce our input file "r1cadd.txt":
# We used the median score per codon; that is 3 scores per nucleotide and 3 nucleotides per codon.
# This gave median Phred-scaled score per codon / residue.
# E.g. 
#Chrom	Pos	PHRED1	PHRED2	PHRED3
#11	36594855	22.3	18.81	22.4
#11	36594856	25.3	23.6	24.6
#11	36594857	24.8	24.3	24.5
#Median PHRED=24.3

# We took the median over a 3 nucleotide window and then only took every 3rd line starting at nucleotide 1 to produce the input data with the correct reading frame:
#awk 'NR == 1 || NR % 3 == 1' r1cadd_median.txt > r1cadd.amino.txt
# r1cadd.txt is output which lists the median PHRED score per codon.

# Dataset: The median PHRED-scaled CADD score per amino acid in RAG1 (generated median as explained in Methods "data processing".
r1cadd <- read.csv(file = "./RAG1.cadd.amino.csv",sep=",",header=T)
mrf.cadd.hgmd <- r1cadd %>% full_join(joined.dedup.short)
mrf.cadd.hgmd1 <- na.omit(mrf.cadd.hgmd, c(hmbd))

# We want to count the scores of MRF and CADD for variants, both reported and unreported, on hgmd. Since the number of reported variants is low compared to the total possible, we can't visually compare the count. Instead the stat density is used to make it easy to view and logically compare.

p6 <-
  mrf.cadd.hgmd1 %>% 
  ggplot(aes(x=Median_cadd, stat(density),  fill=hgmd))+
  geom_density(alpha=0.5, stat = "bin" )+ 
  theme(legend.position="none", panel.background = element_rect("#F7F7F7"))+
  labs(x = "Phred CADD")+
  scale_fill_manual(values=c("#D61A46", "#E1D4F7"))+
  ggtitle("(i)")

p7 <-
  mrf.cadd.hgmd1 %>% 
  ggplot(aes(x=Mutation_rate_residue_frequency, stat(density), fill=hgmd))+
  geom_density(alpha=0.5, stat = "bin" )+ 
  theme(legend.position="none", panel.background = element_rect("#F7F7F7"))+
  labs(x = "MRF score")+
  scale_fill_manual(values=c("#D61A46", "#E1D4F7"), 
                    name = "Reported pathogenic variant",
                    labels = c("Pathogenic", "Unreported"))+
  ggtitle("(ii)")
#grid.arrange(p6, p7, ncol=1)

#Question - 
#  We first produced a list of top candidates and were then asked what % of those are now found in the hgmd database?

# Count mrf categories
mrf.hgmd.count  <- mrf.cadd.hgmd1 %>% 
  dplyr::group_by(Mutation_rate_residue_frequency, Gene, hgmd) %>%
  dplyr::summarise(count = length(Mutation_rate_residue_frequency))

# Round the decimal places to 4 digits for graphing on x-axis
mrf.hgmd.count[, "mrf_rounded"] <- 
  round(mrf.hgmd.count[,"Mutation_rate_residue_frequency"], digits=4)

# Graph the activity measured per mrf category group. Untested variants are shown in grey.
p7b <-
  ggplot() + 
  geom_bar(data = mrf.hgmd.count %>% 
  # convert the MRF score labels fron values to characters to remove x-axis gaps
             dplyr::mutate(mrf_char = as.character(mrf_rounded)), 
           aes(mrf_char, count, fill = hgmd), 
           width = 1, stat = "identity", position = "stack", alpha=0.5) + 
  # plot the same bars with no fill (alpha=0), but colour black to provide bar boarders
  geom_bar(data = mrf.hgmd.count %>%
             dplyr::mutate(mrf_char = as.character(mrf_rounded)), 
          aes(mrf_char, count, fill = hgmd), 
          width = 1, stat = "identity", position = "stack", colour="black", alpha=0)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect("#F7F7F7"),
        legend.position="right")+
  labs(x ="MRF score",
       y = "Number of residues")+
    #annotate("text", x = c(2,10), y = c(20,25), label = c("label 1", "label 2"))+
  annotate("text", x = c(10,18.1), y = c(100,70), 
           label = c((expression(decreasing %<->% increasing)), "36%"))+
scale_fill_manual(values=c("#D61A46", "#E1D4F7"), 
                  name = "HGMD reported \nvariant",
                  labels = c("Pathogenic", "Unreported"))+
  ggtitle("(iii)")

grid.arrange(p6, p7, p7b, ncol=1)
#grid.arrange(p6, p7, p7b, layout_matrix = rbind(c(1,3),c(2,3)))
# size 7.1 x 9

# AUC calculated using pdl6 and pdl7; fill, x (for MRF or CADD score), density.
# Percentage AUC overlap difference above intersect >0.0409 and >22.84 for MRF and CADD, respectively.
# p6.ld <- layer_data(p6, 1)
# p7.ld <- layer_data(p7, 1)

################################################
####  Figure E1                             ####
################################################
####  Reillustration of the data used in    ####
####  Figure 1, now showing the 1% average  ####
####  as a line instead of heatmap, useful  ####
####  if one wants to visualised a          ####
####  horizontal threshold cut, e.g. >0.02  ####
################################################

ggplot(mrf, aes(x=Residue, y=Average_over_1pc, color=Average_over_1pc))+
  geom_line()+
  # split into two separate panels
  facet_wrap(~Gene, nrow =2, scales="free")+ 
  # label the axis and legend
  labs(x = "Residue number",
       y = "MRF",
       color = "MRF with \n1% average \n")+
  # set the colors used  
  scale_color_gradient(low="white", high="red2")+
  # turn down the background darkness and set the legend font size
  theme(panel.background = element_rect("#F7F7F7"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=9))

################################################
####  Figure E2                             ####
################################################
####  "MRF likelihood score versus known    ####
####  functional activity."                 ####
################################################

# We show here all variants that we know to be tested to date using functional assays.
# The bottom panels show that most mutations tested show severe loss of protein function, while the top panels shown how likely each mutation is predicted to occur in humans. 

# Dataset: Prediction MRF and measured recombination by Lee, Tirosh, and Lawless.
mrf.v.activity <- read.csv(file = "./Predicted_measured.csv",sep=",",header=T)
# pc_acuracy <- read.csv(file = "./Predicted_measured/pc_acuracy.csv",sep=",",header=T)

# Set the order to rank by mrf score 
mrf.v.activity.rank <- mrf.v.activity[order(mrf.v.activity$Gene,mrf.v.activity$source,mrf.v.activity$value),]

mrf.v.activity.rank$variant <- factor(mrf.v.activity.rank$variant, levels = rev(unique(mrf.v.activity.rank$variant)), ordered=TRUE)

# Rename the "source" label with something meaningful for the facet labels.
levels(mrf.v.activity.rank$source) <- c("Mutation likelihood\n(% max MRF score)", "Loss of recombination activity\n(% maf loss of function)")

ggplot() + 
  geom_bar(data = mrf.v.activity.rank,
           aes(variant, value, fill = value),
           colour = "slateblue4",
           width = 1, stat = "identity", position = "stack") + 
  scale_fill_gradient(low="cornsilk1", high="thistle", limits=c(0, 100), name = "Percentage\nmaximum") + 
  facet_grid(source~Gene, scales = "free")+
  scale_y_continuous(limits = c(0, 100))+ 
  theme(panel.background = element_rect("#F7F7F7"),
        legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x ="Individual variant residue per bar",
       y = "Percentage of maximum score")
# size 5 x 9

################################################
####  End of Analysis                       ####
################################################
####  End of Analysis                       ####
################################################
####  End of Analysis                       ####
################################################
####  End of Analysis                       ####
################################################

