# Clone the current GitHub repository and specify the path to the following supplementary tables
library(openxlsx)
SUPP = list(
  table.1 = "./6. Supplementary/data/Supplementary Table 1 - Participants Clinicopathological Characteristics.xlsx",
  table.3 = "./6. Supplementary/data/Supplementary Table 3 - Endophenotypes Associations.xlsx",
  table.5 = "./6. Supplementary/data/Supplementary Table 5 - BEYOND analysis results.xlsx")



# ----------------------------------------------------- #
#           Subpopulation proportion matrix             # 
# ----------------------------------------------------- #
# Create and load a conda environment with the necessary R and python packages
source("2. Cell-type analysis/utils/load.code.env.R")
library(anndata)
library(dplyr)
library(tibble)

proportions <- read.xlsx(SUPP$table.3, "snRNA-seq proportions") %>% column_to_rownames("individualID")
clinical.information <- read.xlsx(SUPP$table.1, "Discovery cohort") %>% column_to_rownames("individualID")
qcs <- read.xlsx(SUPP$table.3, "Participant inclusion QCs") %>% filter(Final_QC == "Pass") %>% column_to_rownames("individualID")

# We represent the cellular landscape using an AnnData object. In its core `data$X` is the subpopulation proportion matrix.
# This representation provides us a data structure into which we can load participant- (rows) and subpopulation (column) information
data <- AnnData(
  X = proportions,
  layers = list(sqrt.prev = sqrt(proportions)),
  obsm = list(
    QCs = qcs[rownames(proportions), ],
    meta.data = clinical.information[rownames(proportions), ])
)

# The following `cogdx_ad` is derived from the `cogdx` characterization and is used for subpopulation-endophenotype associations (while excluding MCI or AD with another cause of CI)
data$obsm$meta.data$cogdx_ad = as.numeric(recode(data$obsm$meta.data$cogdx, "1"="1","2"="2", "3"=NA_character_, "4"="3", "5"=NA_character_,"6"=NA_character_))

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(proportions, clinical.information, qcs)


# ----------------------------------------------------- #
#      Endophenotype associations & meta-analysis       # 
# ----------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$uns$celmod <- list()
data$uns$celmod$predicted.proportions <- read.xlsx(SUPP$table.3, "CelMod predicted proportions")
data$uns$celmod$avg.predicted.prop <- 
  data$uns$celmod$predicted.proportions %>% py_to_r %>% split(., .$set) %>% 
  lapply(., function(df) dcast(df, individualID~subpopulation, value.var = "sqrt.prev_mean") %>% 
           column_to_rownames("individualID") %>%
           `[`(,colnames(data$X)))
data$uns$celmod$test.corrs    <- read.xlsx(SUPP$table.3, "CelMod correlations") %>% column_to_rownames("subpopulation")
data$uns$celmod$celmod.states <- data$uns$celmod$test.corrs %>% py_to_r %>% filter(adj.pval < .01 & corr > 0) %>% rownames()
data$uns$celmod$shared.donors <- data$uns$celmod$avg.predicted.prop$train$index

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")


data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")
associations <- read.xlsx(SUPP$table.3, "Endophenotype associations") %>% split(., .$cohort)

data$uns$trait.analysis <- list(
  snuc = associations$discovery, 
  celmod = associations$replication)

data$uns$trait.analysis$meta.analysis <- 
  merge(py_to_r(data$uns$trait.analysis$snuc),
        py_to_r(data$uns$trait.analysis$celmod),
        by = c("trait","state"),
        suffixes = c(".sc",".b"),
        all.x = T) %>% 
  merge(., py_to_r(data$uns$celmod$test.corrs) %>% `colnames<-`(paste0(colnames(.),".celmod")),
        by.x = "state",
        by.y = "row.names") %>% arrange(-corr.celmod) %>% 
  merge(., read.xlsx(SUPP$table.3, "Endophenotype meta-analysis") %>% rename("state"="subpopulation") %>% select(-ends_with(".sc"), -ends_with(".b")), by=c("trait","state"))

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(associations)


# ----------------------------------------------------- #
#         BEYOND - Cellular landscape embedding         # 
# ----------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$obs$clusters <- read.xlsx(SUPP$table.5, "3D Landscape embedding") %>% 
  column_to_rownames("individualID") %>%
  `[`(rownames(data), "cluster")
data$obs$core <- !data$obs$clusters %in% c(9,10)

data$obsm$X_all_3d_phate <- read.xlsx(SUPP$table.5, "3D Landscape embedding") %>% 
  column_to_rownames("individualID") %>% 
  select(-cluster) %>% `[`(rownames(data), )
data$obsm$X_core_phate <- read.xlsx(SUPP$table.5, "2D Landscape embedding") %>% 
  column_to_rownames("individualID") %>% 
  select(-cluster) %>% `[`(rownames(data), )

# Compute embedding local density for smoothened landscape plots
source("4. BEYOND/utils/utils.R")
data$obsp <- list()
for(e in c("X_all_3d_phate","X_core_phate")) 
  data$obsp[[paste0("similarity_", e)]] <- embedding.similarity(data$obsm[[e]], knn = 5)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")


# ----------------------------------------------------- #
#       BEYOND - Cellular trajectories & dynamics       # 
# ----------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

palantir <- read.xlsx(SUPP$table.5, "Palantir Trajectories") %>% column_to_rownames("individualID") %>% `[`(rownames(data), )
via <- read.xlsx(SUPP$table.5, "VIA Trajectories") %>% column_to_rownames("individualID") %>% `[`(rownames(data), )

data$uns$trajectories <- list(
  palantir  = list(
    pseudotime = setNames(palantir$pseudotime, rownames(palantir)),
    branch.probs = palantir[, c("prAD", "ABA")],
    terminals = palantir %>% filter(!is.na(terminal)) %>% 
      select(pseudotime, traj=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("traj"),
    user.root = palantir[!is.na(palantir$root), "root"]
  ),
  
  via = list(
    pseudotime = setNames(via$pseudotime, rownames(via)),
    branch.probs = via[, c("prAD.like", "ABA.like", "Trajectory.3")] %>% `colnames<-`(gsub(".like", "", colnames(.))),
    terminals = via %>% filter(!is.na(terminal)) %>% 
      select(pseudotime, traj=terminal) %>% rownames_to_column("terminal") %>% column_to_rownames("traj"),
    user.root = via[!is.na(via$root), "root"]
  )
)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(palantir, via)


data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

data$uns$trajectories$palantir$dynamics = list(
  fitted.vals = read.xlsx(SUPP$table.5, "Dynamics fitted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(y=value, x=pseudotime, X.weights.=weight, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory),
  
  pred.vals = read.xlsx(SUPP$table.5, "Dynamics predicted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(x=pseudotime, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory)
)

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")


# ----------------------------------------------------- #
#                 BEYOND - Communities                  # 
# ----------------------------------------------------- #
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

# Append sub-population to community partition
data$var <- read.xlsx(SUPP$table.5, "Subpopulation communities") %>% 
  column_to_rownames("subpopulation") %>% `[`(colnames(data), )

# Append participants' community proportions
data$obsm$communities <- read.xlsx(SUPP$table.5, "Participant comm. prop.") %>% 
  column_to_rownames("individualID") %>% `[`(rownames(data), c("C1", "C2", "C3"))
data$obsm$sub.communities <- read.xlsx(SUPP$table.5, "Participant comm. prop.") %>% 
  column_to_rownames("individualID") %>% `[`(rownames(data),) %>% select(-C1, -C2, -C3)

# Append community dynamics
data$uns$communities$dynamics = list(
  fitted.vals = read.xlsx(SUPP$table.5, "Dynamics fitted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(y=value, x=pseudotime, X.weights.=weight, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory),
  
  pred.vals = read.xlsx(SUPP$table.5, "Dynamics predicted values") %>% 
    filter(!feature %in% c("C1", "C2", "C3", "C4", "C5", "C1.1", "C1.2", "C2.3", "C2.2", "C2.1")) %>% 
    mutate(fit_sd=NA, se.fit_sd=NA) %>%
    select(x=pseudotime, fit, se.fit, fit_sd, se.fit_sd, feature, trajectory)
)

# Append endophenotype associations
data$uns$communities$trait.association <- read.xlsx(SUPP$table.5, "Community endophe. assoc.")

anndata::write_h5ad(data, "2. Cell-type analysis/data/subpopulation.proportions.h5ad")
rm(data)


# ----------------------------------------------------- #
#             Plotting manuscript figures               # 
# ----------------------------------------------------- #
source("5. Manuscript code/utils.R")
data <- anndata::read_h5ad("2. Cell-type analysis/data/subpopulation.proportions.h5ad")

landscape <- plot_grid(
  ggdraw() + draw_label("Trajectories in Cellular Landscape"),
  plot.landscape(data$uns$trajectories$palantir$branch.probs %>% # 手动加载5. Manuscript code/utils.R中Plotting Utils Functions
                   py_to_r %>% mutate(diff=prAD-ABA) %>% dplyr::select(diff),
                 smoothened = F,
                 size = 2,
                 cols = green2purple.less.white,
                 legend.position = "right", show.feature.name = FALSE),
  ncol=1,
  rel_heights = c(.1, 1))

dynamics <- plot_grid(
  ggdraw() + draw_label("Subpopulation Dynamics Along Trajectories"),
  plot_grid(
    plot.dynamics(features = c("Oli.7","Ast.10"), 
                  cols = list(Oli.7="olivedrab4", Ast.10="darkorchid4"), 
                  label = T, legend.position = "none", include.points = TRUE, overlap.pseudotime = .1) +
      labs(x="Pseudotime", y="Proportion") + 
      theme(strip.background = element_rect(fill="lightgrey")),
    plot.dynamics(features = c("Ast.5","OPC.3"), 
                  cols = list(Ast.5="darkorchid4",OPC.3="springgreen4"), 
                  label = T, legend.position = "none", include.points = TRUE, overlap.pseudotime = .1) +
      labs(x="Pseudotime", y="Proportion") + 
      theme(strip.background = element_rect(fill="lightgrey")),
    ncol=1),
  ncol=1,
  rel_heights = c(.1, 1))

associations <- plot_grid(
  ggdraw() + draw_label("Endophenotype Associations"),
  plot.trait.associations(py_to_r(data$uns$trait.analysis$snuc), 
                          params = names(AD.traits),
                          column_title="Discovery (snRNA-seq)",
                          column_labels = c("A", "T", "C"),
                          column_names_rot = 0,
                          column_names_centered = T,
                          row_names_side = "left",
                          use_raster=T,
                          border=T,
                          raster_quality = 10) %>% draw %>% grid.grabExpr(),
  ncol=1,
  rel_heights = c(.1, 1))

landscape + dynamics + associations


