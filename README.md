---
title: "Generalized Linear Model with edgeR"
author: "Jiacheng Ding"
date: "2025-04-10"
output:
  html_document:
  toc: true
toc_depth: 2
toc_float: true
number_sections: true
theme: readable
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
```

```{r load-packages}
library(edgeR)
library(ggplot2)
library(compcodeR)
library(dplyr)
library(tidyverse)
library(MASS)
set.seed(123456)
```

## Load bulk RNA-seq Data

```{r load-data}
# setwd removed for portability. Make sure file path is correct.
if (!file.exists('combined_counts.txt')) stop('File not found: combined_counts.txt')
rna_data <- read.table('combined_counts.txt', header = TRUE, check.names = FALSE) %>%
  tibble::column_to_rownames("Geneid")
head(rna_data)
```

## Poisson vs Negative Binomial

```{r poisson-nb-comparison}
# ==== Negative binomal vs Poisson ====
real_counts <- rna_data + 1 # add 1 to avoid log infinite issue
# get per gene mean expression
real_mu <- rowMeans(real_counts)
# get per gene variance
real_var <- apply(real_counts, 1, var)

# Simulate Poisson and NB
n_genes <- nrow(real_counts)
n_samples <- ncol(real_counts)

# Use real means for poisson simulation
mu_pois <- matrix(rep(real_mu, n_samples), nrow = n_genes)
pois_counts <- matrix(rpois(n_genes * n_samples, lambda = mu_pois), nrow = n_genes)

# simplified NB simulation with fixed dispersion across genes
# Use real means for NB simulation
dispersion <- 0.4
nb_counts <- matrix(rnbinom(n_genes * n_samples, 
                            mu = real_mu, 
                            size = 1/dispersion), nrow = n_genes)

# Compute variance for simulated poisson & NB rna counts
pois_var <- apply(pois_counts, 1, var)
nb_var <- apply(nb_counts, 1, var)

# visualize mean expression vs variance at log scale
df_var <- data.frame(
  mean = c(real_mu, real_mu, real_mu),
  variance = c(real_var + 1, pois_var + 1, nb_var + 1), # add 1 to avoid log infinite issue
  model = rep(c("Real Data", "Poisson", "NegBinomial"), each = n_genes)
)

ggplot(df_var, aes(x = mean, y = variance, color = model)) +
  geom_point(alpha = 0.1, size = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Mean Expression (log10)",
    y = "Variance (log10)",
    color = "Model"
  )

```

## Dispersion and Variance

```{r dispersion-variance}
# ==== dispersion & variance ====
# in NB dispersion and variance are linked, however different concepts
# mathematically, for NB, var = u + u^2 * disp, where u is mean
# variance (how much read counts varies from mean) is high at high counts
# However dispersion (driver of var) is high at low counts

# Here we simulate some rna seq data to build up ground truth, 
# using library(compcodeR)
# https://www.rdocumentation.org/packages/compcodeR/versions/1.8.2/topics/generateSyntheticData
sim_DE <- generateSyntheticData(
  dataset = "Validation data",
  n.vars = 10000, # assuming 10,000 genes
  samples.per.cond = 4, # 4 samples per group
  n.diffexp = 0,  # 5% genes are truly differentiated
  fraction.upregulated = 0.5, # 50% genes are upregulated in control group
  dispersions = 'auto', # use pre-trained rna seq dispersion
  effect.size = 4, # log2fold change 2
  repl.id = 1 # simulation id
)

# similarly, we inspect if simulated rna data follows NB
sim_counts <- sim_DE@count.matrix + 1 # add 1 to avoid log infinite issue
# get per gene mean expression
sim_mu <- rowMeans(sim_counts)
# get per gene variance
sim_var <- apply(sim_counts, 1, var)

# Simulate Poisson and NB
n_genes_sim <- nrow(sim_counts)
n_samples_sim <- ncol(sim_counts)

# Use sim means for poisson simulation
mu_pois_4sim <- matrix(rep(sim_mu, n_samples_sim), nrow = n_genes_sim)
pois_counts_4sim <- matrix(rpois(n_genes_sim * n_samples_sim, lambda = mu_pois_4sim), nrow = n_genes_sim)

# simplified NB simulation with fixed dispersion across genes
# Use sim means for NB simulation
dispersion <- 0.4
nb_counts4sim <- matrix(rnbinom(n_genes_sim * n_samples_sim, 
                            mu = sim_mu, 
                            size = 1/dispersion), nrow = n_genes_sim)

# Compute variance for simulated poisson & NB rna counts
pois_var_sim <- apply(pois_counts_4sim, 1, var)
nb_var_sim <- apply(nb_counts4sim, 1, var)

# visualize mean expression vs variance at log scale
df_var_sim <- data.frame(
  mean = c(sim_mu, sim_mu, sim_mu),
  variance = c(sim_var+1, pois_var_sim+1, nb_var_sim+1), # add 1 to avoid log infinite
  model = rep(c("Simulated Data", "Poisson", "NegBinomial"), each = n_genes_sim)
)

ggplot(df_var_sim, aes(x = mean, y = variance, color = model)) +
  geom_point(alpha = 0.1, size = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Mean Expression (log10)",
    y = "Variance (log10)",
    color = "Model"
  )

# we define a function to conveniently extract information from edgeR estimation
# e.g. dispersion, fitted value, residual, etc.
prepare_dge <- function(sim, label, norm_method) {
  group <- factor(sim@sample.annotations$condition) # simulation grouping
  design <- model.matrix(~group) # design matrix
  y <- DGEList(counts = sim@count.matrix, group = group, design = design)
  y$sample <- colnames(sim@count.matrix)
  y$group <- group
  y$dataset <- label
  y$true_lib_size <- colSums(y$counts)
  y <- calcNormFactors(y, method = norm_method)
  y <- estimateDisp(y, tagwise = T) # estimate tag wise dispersion
  return(y)
}

dege_sim <- prepare_dge(sim_DE, 'simulation','TMM')

# How to understand BCV or biological coefficient of variation
# We know CV is SE / u, while SE = sqrt(Var); Var = u + u^2*disp
# Hence CV^2 = 1/u + disp
# 1/u is the noise (1/Var) of poisson model: sequencing sampling variability 
# excess variance (disp) accounts for biological variance 
# BCV = sqrt(disp)
# more to read edgeR tutorial page 18
# based on the above relationship: BCV = sqrt(disp)
# we can generate BCV plot by manually pulling out dispersion and avg logcpm
# so bcv is essentially dispersion
# you can read more about likelihood approach for dispersion estimation
# edgeR tutorial - page 20 - 22
avg_logcpm <- dege_sim$AveLogCPM
bcv <- sqrt(dege_sim$tagwise.dispersion)
par(mfrow = c(1, 2))
# plot BCV
plotBCV(dege_sim)
# plot manually derived bcv
plot(x = avg_logcpm, y = bcv, xlab = 'Avg logCPM', ylab = 'BCV',cex =0.1)
par(mfrow = c(1, 1))

# For the purpose of better realizing the difference between disp and variance
# we visualize dispersion and variance together
# again Var = u + u^2 * disp
# we can also directly inspect relationship between dispersion and avg logcpm
disp <- dege_sim$tagwise.dispersion
par(mfrow = c(1, 2))
plot(x = avg_logcpm, y = disp, xlab = 'Avg logCPM', ylab = 'Dispersion',cex =0.1)
mu <- 2^avg_logcpm
disp <- dege_sim$tagwise.dispersion
variance <- mu + disp * mu^2
plot(avg_logcpm, log2(variance + 1),
     xlab = "Avg logCPM",
     ylab = "log2(Variance)", cex = 0.1)
abline(0, 1, col = "black", lty = 2)
par(mfrow = c(1, 1))

# we can inspect how edgeR estimate dispersion compared to observed disp and true disp
# observed dispersion
disp_est <- (sim_var - sim_mu) / (sim_mu^2)
disp_true <- sim_DE@variable.annotations$truedispersions.S1
disp <- dege_sim$tagwise.dispersion
par(mfrow = c(1, 3))
plot(disp_est, disp,
     xlab = "Observed dispersion",
     ylab = "edgeR stabilized dispersion", xlim = c(0,10), cex = 0.1)
abline(0, 1, col = "black", lty = 2)
plot(disp_true, disp,
     xlab = "True dispersion",
     ylab = "edgeR stabilized dispersion", xlim = c(0,10),cex = 0.1)
abline(0, 1, col = "black", lty = 2)
plot(disp_true, disp_est,
     xlab = "True dispersion",
     ylab = "Observed dispersion", xlim = c(0,10), cex = 0.1)
abline(0, 1, col = "black", lty = 2)
par(mfrow = c(1, 1))
```

## Normalization Importance

```{r normalization-check}
# ==== Normalization ====
# Here we simulate rna data with 0 differential genes 
# Intentionally set 10 fold difference between group 1 and group 2.
sim_noDE <- generateSyntheticData(
  dataset = "noDE",
  n.vars = 10000,
  samples.per.cond = 4,
  n.diffexp = 0,
  fraction.upregulated = 0.5,
  dispersions = 'auto',
  effect.size = 2,
  repl.id = 2,
  seqdepth = c(rep(1e7,4), rep(1e8,4)),
  filter.threshold.total = 0
)
# Run de test w/o normalization
dge_noDE_broken <- DGEList(counts = sim_noDE@count.matrix, group = factor(sim_noDE@sample.annotations$condition))
dge_noDE_broken$group <- factor(sim_noDE@sample.annotations$condition)
# Force effective lib size to be same = libSize * normFactor
dge_noDE_broken$samples$lib.size <- rep(1, 8)  # fake lib.size
dge_noDE_broken$samples$norm.factors <- c(1, 8)
dge_noDE_broken <- estimateDisp(dge_noDE_broken)
dge_noDE_broken_ex <- exactTest(dge_noDE_broken) %>% topTags(n = nrow(dge_noDE_broken$counts))

fp_p_no_norm <- dge_noDE_broken_ex$table %>% dplyr::filter(PValue < 0.05) %>% nrow()
fp_p_no_norm <- round(fp_p_no_norm, 4)
fp_fdr_no_norm <- dge_noDE_broken_ex$table %>% dplyr::filter(FDR < 0.05) %>% nrow()
round(fp_fdr_no_norm, 4)

# check false positive discovery
print(paste0('DE genes (no normalization & p-val 5%) detected:', fp_p_no_norm))
print(paste0('DE genes (no normalization & FDR 5%) detected:',fp_fdr_no_norm))

# run regular de test with TMM normalization
dge_noDE_TMM <- prepare_dge(sim_noDE, "low-DE", 'TMM')
dge_noDE_TMM$samples$norm.factors
dge_noDE_TMM_ex <- exactTest(dge_noDE_TMM) %>% topTags(n = nrow(dge_noDE_TMM$counts))

fp_p_tmm <- dge_noDE_TMM_ex$table %>% dplyr::filter(PValue < 0.05) %>% nrow()
fp_p_tmm <- round(fp_p_tmm, 4)
fp_fdr_tmm <- dge_noDE_TMM_ex$table %>% dplyr::filter(FDR < 0.05) %>% nrow()
fp_fdr_tmm <- round(fp_fdr_tmm, 4)

# check false positive discovery 
print(paste0('DE genes (TMM normalization; p-val 5%) detected:',fp_p_tmm))
print(paste0('DE genes (TMM normalization; FDR 5%) detected:',fp_fdr_tmm))

# check MDS plot (notice % variation explained)
plotMDS(dge_noDE_broken, col = as.numeric(dge_noDE_broken$group),
        main = "Before normalization",top=2000, dim.plot = c(1,2), gene.selection = 'common')
plotMDS(dge_noDE_TMM, col = as.numeric(dge_noDE_TMM$group),
        main = "After TMM normalization",top=2000, dim.plot = c(1,2), gene.selection = 'common')

# check p value inflation
# p-value should be uniformly distributed under null H
hist(dge_noDE_broken_ex$table$PValue, breaks = 50, main = "P-value distribution without normalization")
hist(dge_noDE_TMM_ex$table$PValue, breaks = 50, main = "P-value distribution with TMM")
```

## Linear Model Fitting in edgeR

```{r glm-fit}
# ==== GLM vs CLM for RNA data ====
# CLM assumes normal distribution and assumes constant variance,
# in which we have known as false for rna seq
# for ran seq data we need account for over-dispersion (mean-variance relationship)
# use GLM with NB frame to model mean gene expression and dispersion

# Here we compare CLM and GLM under the same NB assumption
n_genes_glm <- 100
x <- seq(0, 5, length.out = n_genes_glm)  # fake effect size
mu <- exp(x)                              # log-linear mean B = log(u) / exp(B) = u
counts_nb <- rnbinom(n_genes_glm, mu = mu, size = 1/0.1)  # NB-distributed counts

df <- data.frame(x = x, y_nb = counts_nb)
fit_glm <- glm.nb(y_nb ~ x, data = df)
fit_clm <- lm(y_nb ~ x, data = df)

df$fit_glm <- predict(fit_glm, type = "response")
df$fit_clm <- predict(fit_clm)

# visualize CLM and GLM
# agian CLM doesn't account for overdispersion, which assumes fixed variance
ggplot(df, aes(x = x)) +
  geom_point(aes(y = y_nb), color = "black", size = 2, alpha = 0.6) +
  geom_line(aes(y = fit_glm), color = "red", size = 1.2) +
  geom_line(aes(y = fit_clm), color = "blue", linetype = "dashed", size = 1.2) +
  labs(title = "Blue: CLM; Red: GLM",
       y = "Simulated RNA-seq Counts (NB)",
       x = "Effect Size (Beta)") +
  theme_classic()
```

## Design Matrix with/without Intercept

```{r design-matrix}
# ==== design matrix w/wo intercept (common baseline) ====
sim_counts <- sim_DE@count.matrix
sim_info <- sim_DE@sample.annotations
group <- factor(sim_info$condition)

# assuming group 1 Beta 0.1 and group 2 Beta 0.2
# witout intercept group 2 regression coefficient is 0.2
model.matrix(~0 + group)
model.matrix(~0 + group) %*% c(0.1, 0.2)
# with intercept group 2 regression coefficient is 0.2 + 0.1 (baseline)
model.matrix(~group)
model.matrix(~group) %*% c(0.1, 0.2)
```

## Simulating and Correcting Batch Effects

```{r simulate-batch-effects}
# ==== Batch effects ====
# in this section we will explore the effects of simluated batch effects on de discovery
# ---- _simulate no batch effects ----
# Simulate 10 conditions, 10 batches, 20 samples total
# no batch introduced
sim_nobatch <- generateSyntheticData(
  dataset = "batch_example",
  n.vars = 10000,
  samples.per.cond = 10,
  n.diffexp = 500, # 5% de genes
  dispersions = "auto",
  fraction.upregulated = 0.5,
  effect.size = 3,
  seqdepth = 1e7,
  filter.threshold.total = 0,
  between.group.diffdisp = F
)

# check how many true positive genes we can detected when there is no batch effects
sample_info_0batch <- sim_nobatch@sample.annotations
group_0batch <- factor(sample_info_0batch$condition)
design_0batch <- model.matrix(~group_0batch)
dge_0batch <- DGEList(counts = sim_nobatch@count.matrix)
dge_0batch <- calcNormFactors(dge_0batch, 'TMM')
dge_0batch <- estimateDisp(dge_0batch, design = design_0batch, tagwise = T)
fit_0batch <- glmFit(dge_0batch, design_0batch)
fit_test_0batch <- glmLRT(fit_0batch, coef = 2)

# inspect true de gene coefficients
true_DE_0batch <- sim_nobatch@variable.annotations$differential.expression == 1
hist(fit_0batch$coefficients[,2][true_DE_0batch], 
     breaks = 100,
     xlab = 'coefficient', ylab = 'frequency', main = 'true DE gene coefficient')

# ground truth de genes
all_genes <- rownames(sim_nobatch@count.matrix)
true_DE_0batch <- sim_nobatch@variable.annotations$differential.expression
true_DE_genes_0batch <- rownames(sim_nobatch@count.matrix)[true_DE_0batch == 1]
true_nonDE_genes_0batch <- rownames(sim_nobatch@count.matrix)[true_DE_0batch == 0]

# exrtact edgeR de test result
res_0batch <- 
  topTags(fit_test_0batch, 
          n = nrow(fit_test_0batch$table))$table

# use FDR 5%
called_DE_0batch <- rownames(res_0batch)[res_0batch$FDR < 0.05]

# count false positive, true negative, false negative, and true positive
FP_0batch <- sum(called_DE_0batch %in% true_nonDE_genes_0batch)
TN_0batch <- sum(setdiff(all_genes, called_DE_0batch) %in% true_nonDE_genes_0batch)

FN_0batch <- sum(setdiff(all_genes, called_DE_0batch) %in% true_DE_genes_0batch)
TP_0batch <- sum(called_DE_0batch %in% true_DE_genes_0batch)

data.frame(false_positive = FP_0batch, false_negative = FN_0batch, true_negative = TN_0batch, true_positive = TP_0batch)
print(paste0('fasle positive rate: ', round(FP_0batch / (FP_0batch + TN_0batch), 4)))
print(paste0('false negative rate: ', round(FN_0batch / (FN_0batch + TP_0batch), 4)))

# ---- _simulate batch effects ----
# introduce batch effects in two levels
# Batch 1: samples 1–5 (condition A) & Sample 11-15 (condition B)
# Batch 2: samples 6–10 (condition A) & Sample 16-20 (condition B)
# level 1 batch effects: Batch 2 with 2 fold seq depth of Batch 1 (overall expression profile shift)
sim_batch <- generateSyntheticData(
  dataset = "batch_example",
  n.vars = 10000,
  samples.per.cond = 10,
  n.diffexp = 500, # 5% genes consistent with sim_nobatch
  dispersions = "auto",
  fraction.upregulated = 0.5,
  effect.size = 3, # same effect size as sim_nobatch
  seqdepth = c(rep(1e7, 5), rep(1.5e7, 5),rep(1e7, 5), rep(1.5e7, 5)), # introduce level 1 batch effects
  filter.threshold.total = 0,
  between.group.diffdisp = F
)

# level 2 batch effects: composition change
# library size shift isn't enough to mimic real life batch effects
# batch effects introduced by technical reasons might introduce composition changes
# here we spik in a fixed expression value difference on random 1000 (1%) genes 
# e.g. some genes are more efficiently captured by batch 2
# Apply a fold-change shift to simulate a batch-specific bias
sample_info_batch <- sim_batch@sample.annotations
sample_info_batch$batch <- c(rep('batch1', 5), 
                             rep('batch2', 5),
                             rep('batch1', 5),
                             rep('batch2', 5))
biased_genes_1 <- sample(1:500, 500)
biased_genes_2 <- sample(501:10000, 500)
biased_genes <- c(biased_genes_1, biased_genes_2)
batch2_samples <- which(sample_info_batch$batch == "batch2")

# Increase expression consistently
sim_batch@count.matrix[biased_genes, batch2_samples] <-
  sim_batch@count.matrix[biased_genes, batch2_samples] * 4 # introduce noise larger than true effects

dge_batch <- DGEList(counts = sim_batch@count.matrix)
dge_batch$samples$batch <- sample_info_batch$batch
dge_batch$samples$group <- sample_info_batch$condition

# overall library size shift introduced by batch should be corrected by normalization
dge_batch <- calcNormFactors(dge_batch, 'TMM')

# check norm.factors which accounts for effective lib size
# effective lib size = norm.factors * libSize
nf_df <- dge_batch$samples
ggplot(nf_df, aes(x = batch, y = norm.factors, fill = batch)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, color = "black") +
  labs(
    x = "Batch",
    y = "Normalization Factor"
  ) +
  theme_classic() +
  theme(legend.position = "none")

# Design matrices: with and without batch
group_batch <- factor(sample_info_batch$condition)
batch_batch <- factor(sample_info_batch$batch)

design_w_batch <- model.matrix(~group_batch + batch_batch)
design_no_batch <- model.matrix(~group_batch)
# check if matrix full ranked
qr(design_w_batch)$rank == ncol(design_w_batch)

# Estimate dispersion
dge_batch_modeled <- estimateDisp(dge_batch, design = design_w_batch, tagwise = T)
dge_batch_not_modeled <- estimateDisp(dge_batch, design = design_no_batch, tagwise = T)

# check estmiated dispersion difference
disp_df <- data.frame(
  dispersion = c(dge_batch_modeled$tagwise.dispersion,
                 dge_batch_not_modeled$tagwise.dispersion),
  model = rep(c("With batch", "Without batch"),
              each = length(dge_batch_modeled$tagwise.dispersion))
)
disp_df$log_dispersion <- log10(disp_df$dispersion + 1e-6)  # small offset to avoid log(0)

ggplot(disp_df, aes(x = log_dispersion, fill = model)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(
    x = "log10(Dispersion)",
    y = "Number of Genes"
  ) +
  theme_classic() +
  theme(legend.title = element_blank())

# Fit GLMs
fit_batch_modeled <- glmFit(dge_batch_modeled, design_w_batch)
fit_batch_not_modeled <- glmFit(dge_batch_not_modeled, design_no_batch)

# alternative quasi-like model. Allowing extra dispersion
# fit_batch_modeled <- glmQLFit(dge_batch_modeled, design_w_batch)
# fit_batch_not_modeled <- glmQLFit(dge_batch_not_modeled, design_no_batch)

# check if biased gene has positive coefficients (because we intentionally inflated expression)
hist(fit_batch_modeled$coefficients[,3][biased_genes], 
     breaks = 100,
     xlab = 'coefficient', ylab = 'frequency', main = 'baised gene coefficient')

hist(fit_batch_modeled$coefficients[,3][setdiff(1:10000, biased_genes)], 
     breaks = 100,
     xlab = 'coefficient', ylab = 'frequency', main = 'none baised gene coefficient')

# LRT for group coefficient (2nd coefficient: group)
fit_test_batch_modeled <- glmLRT(fit_batch_modeled, coef = 2)
fit_test_batch_not_modeled <- glmLRT(fit_batch_not_modeled, coef = 2)

# alternative QL F-test
# fit_test_batch_modeled <- glmQLFTest(fit_batch_modeled, coef = 2)
# fit_test_batch_not_modeled <- glmQLFTest(fit_batch_not_modeled, coef = 2)

# inspect signal discovery rate
all_genes <- rownames(sim_batch@count.matrix)
true_DE_batch <- sim_batch@variable.annotations$differential.expression
true_DE_genes_batch <- rownames(sim_batch@count.matrix)[true_DE_batch == 1]
true_nonDE_genes_batch <- rownames(sim_batch@count.matrix)[true_DE_batch == 0]

# DE calls from model with and without batch
res_batch_modeled <- topTags(fit_test_batch_modeled, n = nrow(fit_test_batch_modeled$table))$table
res_batch_not_modeled <- topTags(fit_test_batch_not_modeled, n = nrow(fit_test_batch_not_modeled$table))$table

# Define DE calls at FDR < 0.05
called_DE_batch_modeled <- rownames(res_batch_modeled)[res_batch_modeled$FDR < 0.05]
called_DE_batch_not_modeled <- rownames(res_batch_not_modeled)[res_batch_not_modeled$FDR < 0.05]

# With batch
FP_batch_modeled <- sum(called_DE_batch_modeled %in% true_nonDE_genes_batch)
FN_batch_modeled <- sum(setdiff(all_genes, called_DE_batch_modeled) %in% true_nonDE_genes_batch)
TN_batch_modeled <- sum(setdiff(all_genes, called_DE_batch_modeled) %in% true_DE_genes_batch)
TP_batch_modeled <- sum(called_DE_batch_modeled %in% true_DE_genes_batch)
# sum(c(FP_batch_modeled, FN_batch_modeled, TN_batch_modeled, TP_batch_modeled))
# FP_batch_modeled / (FP_batch_modeled + TN_batch_modeled)
# FN_batch_modeled / (FN_batch_modeled + TP_batch_modeled)

# Without batch
FP_batch_not_modeled <- sum(called_DE_batch_not_modeled %in% true_nonDE_genes_batch)
FN_batch_not_modeled <- sum(setdiff(all_genes, called_DE_batch_not_modeled) %in% true_nonDE_genes_batch)
TN_batch_not_modeled <- sum(setdiff(all_genes, called_DE_batch_not_modeled) %in% true_DE_genes_batch)
TP_batch_not_modeled <- sum(called_DE_batch_not_modeled %in% true_DE_genes_batch)
# sum(c(FP_batch_not_modeled, FN_batch_not_modeled, TN_batch_not_modeled, TP_batch_not_modeled))

# FP_batch_not_modeled / (FP_batch_not_modeled + TN_batch_not_modeled)
# FN_batch_not_modeled / (FN_batch_not_modeled + TP_batch_not_modeled)

# batch effects modeling result
data.frame(false_positive = c(FP_batch_modeled, FP_batch_not_modeled),
           false_negative = c(FN_batch_modeled, FN_batch_not_modeled),
           true_negative = c(TN_batch_modeled, TN_batch_not_modeled),
           true_positive = c(TP_batch_modeled, TP_batch_not_modeled),
           model = c('~group + batch','~group'))

```

## Reconstructing GLM & Residual Visualization

```{r reconstruct-glm}

# ==== dispersion, regression coefficient, and residual ====
# in this section we will implement a typical edgeR differential test pipeline
# and try to validate some assumptions manually
sim_counts <- sim_nobatch@count.matrix
sim_dpb <- DGEList(counts = sim_counts[rowSums(sim_counts) >10,])
sim_dpb <- calcNormFactors(sim_dpb, method = 'TMM')

group <- factor(sim_nobatch@sample.annotations$condition)
design <- model.matrix(~group)
sim_dpb <- estimateDisp(sim_dpb, design = design, tagwise = T)
fit <- glmFit(sim_dpb, design = design)
fit.test <- glmLRT(fit)
topTags(fit.test)

# let's try to manually derive u using the knowledge we have known
# log(u) = log(libSize * normaFactor) + X*B
# read more on edgeR tutorial page 21
offset <- log((fit$samples$lib.size) * fit$samples$norm.factors)
xbeta <- fit$design %*% t(fit$coefficients) %>% t()
offset_matrix <- matrix(rep(offset, each = nrow(xbeta)), nrow = nrow(xbeta))
log_mu <- xbeta + offset_matrix
mu_manual <- exp(log_mu)
# ideally fit$fitted.values should be equal to mu
# the minor discrepancy here is likely due to prior.count added to observed data
fit$fitted.values[1:3,1:3]
mu_manual[1:3,1:3]
# we can inspect average error use  mean squared error
mean((mu_manual - fit$fitted.values)^2)

# edgeR does not store residual in output, though it was calculated internally
# see bioconductor for getting residuals
# https://support.bioconductor.org/p/57617/
# example: get pearson residuals
yo <- fit$counts # observation about rna counts
mu <- fit$fitted.values # expectation or model fitted values
phi <- fit$dispersion # equal to sim_dpb$tagwise.dispersion
var <- mu*(1+phi*mu) # variance
# notice here yo - mu is obervation - expectation, which is then normalized to sqrt variance
resid.pearson <- (yo-mu) / sqrt(var)

# reconstruct generalized linear model
# Choose a top DE gene to visualize
top_gene <- rownames(topTags(fit.test, n = 1)$table)[1]
# Pull data for selected gene
obs_counts <- yo[top_gene, ]
fitted_counts <- mu[top_gene, ]
resid_values <- resid.pearson[top_gene, ]
groups <- group
plot_df <- data.frame(
  Sample = colnames(sim_counts),
  Observed = obs_counts,
  Fitted = fitted_counts,
  Residual = resid_values,
  Group = groups
)

# for plot purpose
# make group jitter consistent between residual line and dots
plot_df$Group_num <- as.numeric(plot_df$Group)
plot_df$x_jittered <- jitter(plot_df$Group_num, amount = 0.1)

ggplot(plot_df, aes(x = x_jittered)) +
  geom_segment(aes(y = Fitted, yend = Observed, xend = x_jittered), 
               linetype = "dashed", color = "gray") +
  geom_point(aes(y = Observed, color = "Observed"), size = 3) +
  geom_point(aes(y = Fitted, color = "Fitted"), shape = 17, size = 3) +
  scale_color_manual(values = c("Observed" = "green", "Fitted" = "red")) +
  labs(
    title = paste("GLM Fit and Residuals for", top_gene),
    y = "Count",
    color = "Legend"
  ) +
  theme_classic()

# we can also add the regression line the plot using cofficients
coef_gene <- fit$coefficients[top_gene, ]
intercept <- coef_gene[1]
group_coef <- coef_gene[2]

# Calculate log(mu) = X * beta + offset (mean offset for each group)
# Use average offset per group for visualization
offsets <- offset  # vector of log(libSize * norm.factor)
group_labels <- levels(group)

# Compute average offset per group
avg_offset_by_group <- tapply(offsets, group, mean)

# Predict fitted values per group using log(mu) = intercept + group_coef * group + offset
group_logmu <- c(
  intercept + avg_offset_by_group["1"],                # baseline
  intercept + group_coef + avg_offset_by_group["2"]    # group2 effect
)
group_mu <- exp(group_logmu)  # back-transform to count scale

# Create a data frame for line
model_line <- data.frame(
  Group = group_labels,
  Predicted = group_mu
)

ggplot(plot_df, aes(x = x_jittered)) +
  geom_segment(aes(y = Fitted, yend = Observed, xend = x_jittered), 
               linetype = "dashed", color = "gray") +
  geom_point(aes(y = Observed, color = "Observed"), size = 3) +
  geom_point(aes(y = Fitted, color = "Fitted"), shape = 17, size = 3) +
  geom_line(data = model_line, aes(x = as.numeric(Group), y = Predicted, group = 1), 
            color = "blue", linewidth = 1.2) +
  scale_color_manual(values = c("Observed" = "green", "Fitted" = "red")) +
  scale_x_continuous(breaks = c(1, 2), labels = levels(plot_df$Group)) +
  labs(
    title = paste("GLM Fit and Residuals for gene:", top_gene),
    x = "Group",
    y = "Count",
    color = "Legend"
  ) +
  theme_classic()

```
