# set a working directory
setwd("~/Desktop/COVID-19_Vaccination_GM")

# source R settings
source("01_SCRIPTS/Settings.R")

# read pre-processed files back in
meta_new <- readRDS("02_RESULTS/Modified_metadata.rds")
Tax_table_new <- readRDS("02_RESULTS/Modified_tax_table.rds")
Otu_table <- readRDS("02_RESULTS/Modified_otu_table.rds")
ps <- readRDS("02_RESULTS/Phyloseq.rds")

# transform takes abundances (not = to 1 within sample) and makes them relative abundances (=1 within sample)
ps_main_composition <- microbiome::transform(ps, "compositional")


# DESeq analysis - main Figures 2A 
# factor sample timepoint for analysis
sample_data(ps_main_composition)$Timepoint <- as.factor(sample_data(ps_main_composition)$Timepoint)

# subset the phyloseq object per cohort + take species level taxonomy for eventual heatmap plot
ps_subset <- subset_samples(ps_main_composition, Cohort %in% "ICP")
tax_table(ps_subset) <- tax_table(ps_subset)[,c(1, 2 ,7)]
ps_subset_DESeq <- tax_glom(ps_subset, taxrank = 'Species', NArm = F)

# pairwise comparison between two things of choice
ps_subset_DESeq_test <- subset_samples(ps_subset_DESeq, Timepoint %in% c("Pre-Dose", "Acute"))

# filter sparse features, with > 95% zeros
ps_subset_DESeq_test <- prune_taxa(rowSums(otu_table(ps_subset_DESeq_test) == 0) < ncol(otu_table(ps_subset_DESeq_test)) * 0.95, ps_subset_DESeq_test)

#convert otu  to integers 
otu_table(ps_subset_DESeq_test) <- otu_table(ps_subset_DESeq_test)*100000

# convert to DESeq factorising on Timepoint
ps_ds = phyloseq_to_deseq2(ps_subset_DESeq_test, ~Timepoint)

# run the DESeq
ds = DESeq(ps_ds, test = "Wald", fitType="local", sfType = "poscounts")
alpha = 0.05
res = results(ds, alpha=alpha, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:35, ]) 

# see what species are significantly different between the two timepoints
taxa_sig_check = res[(res$padj < alpha), ]

# acquire data for the heatmap
ps.taxa.rel <- transform_sample_counts(ps_subset_DESeq_test, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps_subset_DESeq_test)), ps.taxa.rel.sig)

otu_matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(otu_matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))

# Define the annotation color for columns and rows
annotation_col = data.frame(`Timepoint` = as.factor(metadata_sub$Timepoint), 
                            check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"]))
rownames(annotation_row) = rownames(otu_matrix)

phylum_col = brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)

# Define the phylum names (replace with your actual phylum names)
phylum_names <- c("p__Firmicutes", "p__Bacteroidetes", "p__Proteobacteria",
                  "p__Actinobacteria", "p__Chloroflexi", "p__Verrucomicrobia",
                  "p__Tenericutes", "p__Bacteria phylum incertae sedis", "p__Fusobacteria",
                  "p__Euryarchaeota", "p__Bacteria phylum [Proteobacteria/Firmicutes]", "p__Elusimicrobia",
                  "p__Lentisphaerae", "p__Spirochaetes", "p__Candidatus Melainabacteria",
                  "p__Bacteria phylum [Firmicutes/Actinobacteria]", "p__Synergistetes", "p__Candidatus Saccharibacteria")

# Define the number of phylum
num_phylum <- length(phylum_names)
# Define a color palette using RColorBrewer
color_palette <- getOI(num_phylum)
# Create a named vector with phylum names as names and corresponding colors
phylum_colors <- setNames(color_palette, phylum_names)

ann_colors = list(`Timepoint` = c("Pre-Dose" = "#E69F00", "Acute" = "#CC79A7"),
                  Phylum = phylum_colors)

Fig2A <- pheatmap(log(otu_matrix+1e-05), scale= "none", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors,
                         cluster_rows = F,
                         cluster_cols = T,
                         show_colnames = F,
                         fontsize = 7,
                         heatmap_legend_param = list(title = "log2 Fold Change"))

# pheatmap functions to save a pdf
pdf("03_FIGURES/Figure2A.pdf",  width = 8, height = 6)
print(Fig2A)
dev.off() 

# save otu_table used for DESeq
# saveRDS(otu_table, "00_DATA/Fig2A_data.rds")


# plot fold change on the significantly different species 
plot_sig = cbind(as(taxa_sig_check, "data.frame"), as(tax_table(ps_subset_DESeq_test)[rownames(taxa_sig_check), ], "matrix"))
x = tapply(plot_sig$log2FoldChange, plot_sig$Species, function(x) max(x))
x = sort(x, TRUE)
plot_sig$Species = factor(as.character(plot_sig$Species), levels=names(x))

Fig2B <- ggplot(plot_sig, 
       aes(x=Species, 
           y=log2FoldChange, 
           color=Phylum)) + 
  geom_point(size=6) + 
  theme(axis.text = element_text(size=16), 
             axis.title = element_text(size=16),
             axis.text.x = element_text(angle = 45, hjust=1, size = 16),
             legend.position = "bottom",
             legend.text = element_text(size = 16),
             legend.title = element_text(size=16))

ggsave("03_FIGURES/Figure2B.pdf", Fig2B)

# save otu_table used for DESeq
# saveRDS(plot_sig, "00_DATA/Fig2B_data.rds")


# DESeq analysis - Supp_Figures 2A + 2C = HC cohort

# subset the phyloseq object per cohort + take species level taxonomy for eventual heatmap plot
ps_subset <- subset_samples(ps_main_composition, Cohort %in% "HC")
tax_table(ps_subset) <- tax_table(ps_subset)[,c(1, 2 ,7)]
ps_subset_DESeq <- tax_glom(ps_subset, taxrank = 'Species', NArm = F)

# pairwise comparison between two things of choice
ps_subset_DESeq_test <- subset_samples(ps_subset_DESeq, Timepoint %in% c("Pre-Dose", "Acute"))

# filter sparse features, with > 95% zeros
ps_subset_DESeq_test <- prune_taxa(rowSums(otu_table(ps_subset_DESeq_test) == 0) < ncol(otu_table(ps_subset_DESeq_test)) * 0.95, ps_subset_DESeq_test)

#convert otu  to integers 
otu_table(ps_subset_DESeq_test) <- otu_table(ps_subset_DESeq_test)*100000

# convert to DESeq factorising on Timepoint
ps_ds = phyloseq_to_deseq2(ps_subset_DESeq_test, ~Timepoint)

# run the DESeq
ds = DESeq(ps_ds, test = "Wald", fitType="local", sfType = "poscounts")
alpha = 0.05
res = results(ds, alpha=alpha, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:35, ]) 

# see what species are significantly different between the two timepoints
taxa_sig_check = res[(res$padj < alpha), ]

# acquire data for the heatmap
ps.taxa.rel <- transform_sample_counts(ps_subset_DESeq_test, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps_subset_DESeq_test)), ps.taxa.rel.sig)

otu_matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(otu_matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))

# Define the annotation color for columns and rows
annotation_col = data.frame(`Timepoint` = as.factor(metadata_sub$Timepoint), 
                            check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"]))
rownames(annotation_row) = rownames(otu_matrix)

phylum_col = brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)

# Define the phylum names (replace with your actual phylum names)
phylum_names <- c("p__Firmicutes", "p__Bacteroidetes", "p__Proteobacteria",
                  "p__Actinobacteria", "p__Chloroflexi", "p__Verrucomicrobia",
                  "p__Tenericutes", "p__Bacteria phylum incertae sedis", "p__Fusobacteria",
                  "p__Euryarchaeota", "p__Bacteria phylum [Proteobacteria/Firmicutes]", "p__Elusimicrobia",
                  "p__Lentisphaerae", "p__Spirochaetes", "p__Candidatus Melainabacteria",
                  "p__Bacteria phylum [Firmicutes/Actinobacteria]", "p__Synergistetes", "p__Candidatus Saccharibacteria")

# Define the number of phylum
num_phylum <- length(phylum_names)
# Define a color palette using RColorBrewer
color_palette <- getOI(num_phylum)
# Create a named vector with phylum names as names and corresponding colors
phylum_colors <- setNames(color_palette, phylum_names)

ann_colors = list(`Timepoint` = c("Pre-Dose" = "#E69F00", "Acute" = "#CC79A7"),
                  Phylum = phylum_colors)

Supp_Fig2A <-  pheatmap(log(otu_matrix+1e-05), scale= "none", 
                  annotation_col = annotation_col, 
                  annotation_row = annotation_row, 
                  annotation_colors = ann_colors,
                  cluster_rows = F,
                  cluster_cols = T,
                  show_colnames = F,
                  fontsize = 7,
                  heatmap_legend_param = list(title = "log2 Fold Change"))

# pheatmap functions to save a pdf
pdf("03_FIGURES/Supp_Fig2A.pdf",  width = 8, height = 6)
print(Supp_Fig2A)
dev.off() 

# save otu_table used for DESeq
# saveRDS(otu_table, "00_DATA/Supp_Fig2A_data.rds")

# plot fold change on the significantly different species 
plot_sig = cbind(as(taxa_sig_check, "data.frame"), as(tax_table(ps_subset_DESeq_test)[rownames(taxa_sig_check), ], "matrix"))
x = tapply(plot_sig$log2FoldChange, plot_sig$Species, function(x) max(x))
x = sort(x, TRUE)
plot_sig$Species = factor(as.character(plot_sig$Species), levels=names(x))

Supp_Fig2C <-  ggplot(plot_sig, 
                aes(x=Species, 
                    y=log2FoldChange, 
                    color=Phylum)) + 
  geom_point(size=6) + 
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust=1, size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16))

ggsave("03_FIGURES/Supp_Fig2C.pdf", Supp_Fig2C)

# save otu_table used for DESeq
# saveRDS(plot_sig, "00_DATA/Supp_Fig2C_data.rds")


# DESeq analysis - Supp_Figures 2B + 2D = PID cohort

# subset the phyloseq object per cohort + take species level taxonomy for eventual heatmap plot
ps_subset <- subset_samples(ps_main_composition, Cohort %in% "PID")
tax_table(ps_subset) <- tax_table(ps_subset)[,c(1, 2 ,7)]
ps_subset_DESeq <- tax_glom(ps_subset, taxrank = 'Species', NArm = F)

# pairwise comparison between two things of choice
ps_subset_DESeq_test <- subset_samples(ps_subset_DESeq, Timepoint %in% c("Pre-Dose", "Acute"))

# filter sparse features, with > 95% zeros
ps_subset_DESeq_test <- prune_taxa(rowSums(otu_table(ps_subset_DESeq_test) == 0) < ncol(otu_table(ps_subset_DESeq_test)) * 0.95, ps_subset_DESeq_test)

#convert otu  to integers 
otu_table(ps_subset_DESeq_test) <- otu_table(ps_subset_DESeq_test)*100000

# convert to DESeq factorising on Timepoint
ps_ds = phyloseq_to_deseq2(ps_subset_DESeq_test, ~Timepoint)

# run the DESeq
ds = DESeq(ps_ds, test = "Wald", fitType="local", sfType = "poscounts")
alpha = 0.05
res = results(ds, alpha=alpha, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:35, ]) 

# see what species are significantly different between the two timepoints
taxa_sig_check = res[(res$padj < alpha), ]

# acquire data for the heatmap
ps.taxa.rel <- transform_sample_counts(ps_subset_DESeq_test, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)
ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps_subset_DESeq_test)), ps.taxa.rel.sig)

otu_matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(otu_matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))

# Define the annotation color for columns and rows
annotation_col = data.frame(`Timepoint` = as.factor(metadata_sub$Timepoint), 
                            check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"]))
rownames(annotation_row) = rownames(otu_matrix)

phylum_col = brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)

# Define the phylum names (replace with your actual phylum names)
phylum_names <- c("p__Firmicutes", "p__Bacteroidetes", "p__Proteobacteria",
                  "p__Actinobacteria", "p__Chloroflexi", "p__Verrucomicrobia",
                  "p__Tenericutes", "p__Bacteria phylum incertae sedis", "p__Fusobacteria",
                  "p__Euryarchaeota", "p__Bacteria phylum [Proteobacteria/Firmicutes]", "p__Elusimicrobia",
                  "p__Lentisphaerae", "p__Spirochaetes", "p__Candidatus Melainabacteria",
                  "p__Bacteria phylum [Firmicutes/Actinobacteria]", "p__Synergistetes", "p__Candidatus Saccharibacteria")

# Define the number of phylum
num_phylum <- length(phylum_names)
# Define a color palette using RColorBrewer
color_palette <- getOI(num_phylum)
# Create a named vector with phylum names as names and corresponding colors
phylum_colors <- setNames(color_palette, phylum_names)

ann_colors = list(`Timepoint` = c("Pre-Dose" = "#E69F00", "Acute" = "#CC79A7"),
                  Phylum = phylum_colors)

Supp_Fig2B <- pheatmap(log(otu_matrix+1e-05), scale= "none", 
                        annotation_col = annotation_col, 
                        annotation_row = annotation_row, 
                        annotation_colors = ann_colors,
                        cluster_rows = F,
                        cluster_cols = T,
                        show_colnames = F,
                        fontsize = 7,
                        heatmap_legend_param = list(title = "log2 Fold Change"))

# pheatmap functions to save a pdf
pdf("03_FIGURES/Supp_Fig2B.pdf",  width = 8, height = 6)
print(Supp_Fig2B)
dev.off() 

# save otu_table used for DESeq
# saveRDS(otu_table, "00_DATA/Supp_Fig2B_data.rds")

# plot fold change on the significantly different species 
plot_sig = cbind(as(taxa_sig_check, "data.frame"), as(tax_table(ps_subset_DESeq_test)[rownames(taxa_sig_check), ], "matrix"))
x = tapply(plot_sig$log2FoldChange, plot_sig$Species, function(x) max(x))
x = sort(x, TRUE)
plot_sig$Species = factor(as.character(plot_sig$Species), levels=names(x))

Supp_Fig2D <-  ggplot(plot_sig, 
                      aes(x=Species, 
                          y=log2FoldChange, 
                          color=Phylum)) + 
  geom_point(size=6) + 
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust=1, size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size=16))

ggsave("03_FIGURES/Supp_Fig2D.pdf", Supp_Fig2D)

# save otu_table used for DESeq
# saveRDS(plot_sig, "00_DATA/Supp_Fig2D_data.rds")