# Function to check and install CRAN or GitHub packages
check_install <- function(pkg, from_github = FALSE, repo = NULL, manager = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (from_github) {
      if (manager == "devtools") {
        if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
        devtools::install_github(repo)
      } else {
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
        remotes::install_github(repo)
      }
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Install core management tools
check_install("devtools")
check_install("librarian")
check_install("BiocManager")

# Install GitHub package
check_install("annotables", from_github = TRUE, repo = "stephenturner/annotables")

# Install Bioconductor packages
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}
library(Biobase)

# Exploratory Data Analysis packages
eda_packages <- c("Amelia", "corrplot", "corrgram", "jtools", "ggiraph", 
                  "ggiraphExtra", "arm", "modelsummary")
lapply(eda_packages, check_install)

# Install additional packages used frequently
common_packages <- c("dplyr", "readxl", "pheatmap", "DESeq2", "tidyverse", "ggplot2", 
                     "ggrepel", "ashr", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", 
                     "KEGGREST", "textshape", "svDialogs", "pathview", "ggpubr", "openxlsx")
lapply(common_packages, check_install)

# Load all at once with librarian
librarian::shelf(dplyr, readxl, pheatmap, DESeq2, tidyverse, ggplot2, ggrepel, 
                 annotables, ashr, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, 
                 KEGGREST, textshape, svDialogs, pathview, ggpubr, openxlsx, Amelia, 
                 corrplot, corrgram, jtools, ggiraph, ggiraphExtra, arm, modelsummary, vegan)

setwd("~/Library/CloudStorage/OneDrive-Personal/PD Guelph/pWAT project/proteome mouse")
data_human <- read_excel("pWAT lean obese mouse analysis full.xlsx", sheet = 1)

data.clean <- na.omit(data_human) #clean to 1440
str(data.clean)


### t-test for lean vs obese ###
p.values <- c()
fc.vec <- c()
for (i in seq_len(nrow(data.clean))) {
  t <- t.test(as.numeric(data.clean[i, 2:6]), as.numeric(data.clean[i, 7:11]), var.equal = TRUE) # test groups
  p.values[i] <- t$p.value
  fc <- mean(as.numeric(data.clean[i, 7:11]))/mean(as.numeric(data.clean[i, 2:6])) # calculate fold change 
  fc.vec <- c(fc.vec, fc)
}

sum(p.values < 0.1)
print(p.values)

# adjust p-values for multiple comparison
p.fdr <- p.adjust(p.values, method = "fdr")
sum(p.fdr < 0.1) 

# create dataframe with significant proteins only
p.fdr_01 <- c(p.fdr < 0.1)
clean_stats <- data.clean
clean_stats$fc <- fc.vec
clean_stats$p.fdr <- p.values
significant_proteins <- clean_stats[p.fdr_01,]

#data for volcano plot before putting in the volcano plot
plogvolcano <- -log10(p.fdr)
fcvolcano <- log2(fc.vec)
clean_stats$fcvolc <- fcvolcano
clean_stats$pvolc <- plogvolcano


#Cleaning data for PCA, volcano plot, and heatmap
clean_prot <- clean_stats %>%
  dplyr::select(-ID, -fc, -p.fdr, -fcvolc, -pvolc)


clean_prot <- as.data.frame(lapply(clean_prot, as.numeric), stringsAsFactors = F)
str(clean_prot)

cor(clean_prot)
mean(cor(clean_prot))

#pca
pca <- rda(decostand(clean_prot, method = "hellinger"), scale = TRUE, center = TRUE)
pca

summary(pca)
scores(pca)
pca

#visualizing with biplot
biplot(pca, scaling = "symmetric") #samples arrows with small angles to an axis are highly correlated with that axis

#eigenvalues
pca.data <- as.data.frame(pca$CA$eig)

#results for variables - lipids
pca.data.prot <- as.data.frame(pca$CA$v)
#results for individuals - mouse
pca.data.mouse <- as.data.frame(pca$CA$u)

#including pc1 and pc2 to the table
pca.data.prot$Groups <- c("lean", "lean", "lean", "lean", "lean",
                          "obese", "obese", "obese", "obese", "obese")

head(pca.data.prot)

#plotting PCA
proteomic_pca <- ggplot(pca.data.prot, aes(x=PC1, y=PC2, color=Groups)) +
  geom_point(size=4) + 
  ggtitle("lean vs. obese AT proteome") +
  xlab("PC1 \n (33% variance)") + ylab("PC2 \n (29% variance)") +
  coord_fixed(ratio = 1) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme_bw() + theme(aspect.ratio = 1) + theme(panel.grid = element_blank())
proteomic_pca + 
  geom_text(aes(label = Groups), size = 2)

proteomic_pca + stat_ellipse(geom = "polygon",
                             aes(fill = Groups), 
                             alpha = 0.25)

#done
#now heatmap
install.packages("pheatmap")
library(pheatmap)

colnames(significant_proteins) <- c("ID", 
                                    "lean1", "lean2", "lean3", "lean4", "lean5",
                                    "obese1", "obese2", "obese3", "obese4",
                                    "Accession", "fc", "adj.p")

df_heatmap <- significant_proteins %>%
  dplyr::select(-ID, -Accession, -fc, -adj.p, -14, -15)
df_heatmap <- as.data.frame(lapply(df_heatmap, as.numeric), stringsAsFactors = F)

p <- pheatmap(df_heatmap, scale = "row",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              cluster_cols = F,
              cluster_rows = T)

#heatmap done without names because there are 502 diff expressed proteins

#volcano plot
colnames(clean_stats) <- c("CODE", "ID_HUMAN", "ID", 
                           "lean", "lean", "lean", "lean", "lean",
                           "obese", "obese", "obese", "obese",
                           "Accession", "fc", "adj.p", "fcvolc", "pvolc")


vp <- clean_stats %>%
  dplyr::select(ID, fcvolc, pvolc)

#theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1:1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(0,20,0,0), size = rel(1:1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#okay, now I converted everything, but let's see if it is working
ggplot(data = vp, aes(x = fcvolc, y = pvolc)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = c(0.1), col = "gray", linetype = "dashed") +
  geom_point()


#up and down color coded
vp$diffexpressed <- 'NO'
vp$diffexpressed[vp$fcvolc > 1 & vp$pvolc > 0.05] <- 'UP'
vp$diffexpressed[vp$fcvolc < -1 & vp$pvolc > 0.05] <- 'DOWN'

#pasting the ggplot again
ggplot(data = vp, aes(x = fcvolc, y = pvolc, col = diffexpressed)) + 
  geom_vline(xintercept = c(0), col = 'gray', linetype = 'dashed') +
  #geom_hline(yintercept = c(0.5), col = 'gray', linetype = 'dashed') +
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Non significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 5), xlim = c(-5, 5)) +
  scale_x_continuous(breaks = seq(-10, 10, 1)) +
  labs(color = 'Grouping', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj.p-value"))


#labeling genes in volcano plot
top_upregulated <- vp %>%
  arrange(desc(pvolc)) %>%
  head(400)

#pasting the ggplot again
ggplot(data = vp, aes(x = fcvolc, y = pvolc, col = diffexpressed)) + 
  geom_vline(xintercept = c(0), col = 'gray', linetype = 'dashed') +
  #geom_hline(yintercept = c(0.5), col = 'gray', linetype = 'dashed') +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Non significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(-5, 5)) +
  scale_x_continuous(breaks = seq(-10, 10, 1)) +
  labs(color = 'Grouping', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj.p-value")) +
  geom_text(data = top_upregulated, aes(label = ID), nudge_x = 0.2, nudge_y = 0.2, check_overlap = T)


