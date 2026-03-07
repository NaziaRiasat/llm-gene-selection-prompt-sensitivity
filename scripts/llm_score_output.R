library(jsonlite)
library(dplyr)
library(readr)

#######################################################
###############A) P1:Jaccard vs ground truth (FDR ≤ 0.05)
#######################################################

########STEP 0 — Ground truth (run once)

gt <- read.csv(DESeq2_top200_minimal.csv,
               stringsAsFactors = FALSE)

gt_sig05 <- sort(as.character(gt$gene[!is.na(gt$padj) & gt$padj <= 0.05]))
gt_n05   <- length(gt_sig05)

cat("Ground-truth DE genes (FDR ≤ 0.05):", gt_n05, "\n")

#####STEP 1 — Helper functions (run once)

library(jsonlite)
library(dplyr)

jaccard_set <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

overlap_coef <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  
  length(intersect(a, b)) / min(length(a), length(b))
}

score_llm_folder <- function(folder, gt_sig, gt_n, llm_name) {
  
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 1)
  
  scores <- lapply(files, function(f) {
    out <- fromJSON(f)
    
    llm_genes <- out$significant_genes
    if (is.null(llm_genes)) llm_genes <- character(0)
    if (is.list(llm_genes)) llm_genes <- unlist(llm_genes, use.names=FALSE)
    llm_genes <- sort(unique(as.character(llm_genes)))
    
    n_reported <- suppressWarnings(as.integer(out$n_significant))
    if (is.na(n_reported)) n_reported <- length(llm_genes)
    
    data.frame(
      llm        = llm_name,
      file       = basename(f),
      n_reported = n_reported,
      n_listed   = length(llm_genes),
      n_true     = gt_n,
      precision  = ifelse(length(llm_genes)==0, NA_real_,
                          length(intersect(llm_genes, gt_sig))/length(llm_genes)),
      recall     = ifelse(gt_n==0, NA_real_,
                          length(intersect(llm_genes, gt_sig))/gt_n),
      jaccard    = jaccard_set(llm_genes, gt_sig),
      overlap = overlap_coef(llm_genes, gt_sig),
      exact_set  = setequal(llm_genes, gt_sig),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(scores)
}

######STEP 2a — Run Prompt 1 for each LLM

scores_p1_chatgpt <- score_llm_folder(
  "ChatGPT/prompt_v1",
  gt_sig05, gt_n05,
  "ChatGPT"
)

scores_p1_gemini <- score_llm_folder(
  "Gemini/prompt_v1",
  gt_sig05, gt_n05,
  "Gemini"
)

scores_p1_claude <- score_llm_folder(
  "Claude/prompt_v1",
  gt_sig05, gt_n05,
  "Claude"
)

#######STEP 3 — Combine + summarize (paper-ready)

scores_p1_all <- bind_rows(
  scores_p1_chatgpt,
  scores_p1_gemini,
  scores_p1_claude
)

print(scores_p1_all)

summary_p1 <- scores_p1_all %>%
  group_by(llm) %>%
  summarise(
    runs = n(),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall = mean(recall, na.rm=TRUE),
    mean_jaccard = mean(jaccard, na.rm=TRUE),
    mean_overlap = mean(overlap, na.rm=TRUE),
    exact_match_rate = mean(exact_set)
  )
print(summary_p1)
dir.create("results", showWarnings = FALSE)
#write.csv(scores_p1_all, "results/prompt1_run_level_scores.csv", row.names=FALSE)
#write.csv(summary_p1,   "results/prompt1_summary_by_llm.csv", row.names=FALSE)

###########################################################################
##########B) P2: Jaccard vs ground truth (0.05 < FDR ≤ 0.10)###################
###########################################################################
#Prompt 5 (FDR ≤ 0.10)
gt_sig_p5 <- sort(as.character(gt$gene[!is.na(gt$padj) & gt$padj <= 0.10]))
gt_n_p5   <- length(gt_sig_p5)
cat("Ground-truth genes for Prompt 5 (FDR ≤ 0.10):", gt_n_p5, "\n")

###Run Prompt 5 for all 3 LLMs

scores_p5_chatgpt <- score_llm_folder("ChatGPT/prompt_v5", gt_sig_p5, gt_n_p5, "ChatGPT")
scores_p5_gemini  <- score_llm_folder("Gemini/prompt_v5",  gt_sig_p5, gt_n_p5, "Gemini")
scores_p5_claude  <- score_llm_folder("Claude/prompt_v5",  gt_sig_p5, gt_n_p5, "Claude")

scores_p5_all <- dplyr::bind_rows(scores_p5_chatgpt, scores_p5_gemini, scores_p5_claude)
print(scores_p5_all)

summary_p5 <- scores_p5_all %>%
  group_by(llm) %>%
  summarise(
    runs = n(),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall    = mean(recall, na.rm=TRUE),
    mean_jaccard   = mean(jaccard, na.rm=TRUE),
    exact_match_rate = mean(exact_set)
  )


summary_p5 <- scores_p5_all %>%
  group_by(llm) %>%
  summarise(
    runs = n(),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall = mean(recall, na.rm=TRUE),
    mean_jaccard = mean(jaccard, na.rm=TRUE),
    mean_overlap = mean(overlap, na.rm=TRUE),
    exact_match_rate = mean(exact_set)
  )

print(summary_p5)

#dir.create("results", showWarnings = FALSE)
#write.csv(scores_p5_all, "results/prompt5_run_level_scores.csv", row.names=FALSE)
#write.csv(summary_p5,   "results/prompt5_summary_by_llm.csv", row.names=FALSE)


#######################################
###########C: Promnpt P6#############
#######################################

##########helper functions:

jaccard_set <- function(a, b) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}
overlap_coef <- function(a, b) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / min(length(a), length(b))
}

read_top20 <- function(f, key="top20_genes") {
  out <- fromJSON(f)
  genes <- out[[key]]
  if (is.null(genes)) genes <- character(0)
  if (is.list(genes)) genes <- unlist(genes, use.names = FALSE)
  sort(unique(as.character(genes)))
}

score_folder_top20 <- function(folder, llm_name, gt_top20) {
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 1)
  
  scores <- lapply(files, function(f) {
    llm_genes <- read_top20(f, key="top20_genes")
    
    data.frame(
      llm = llm_name,
      file = basename(f),
      n_listed = length(llm_genes),
      precision = ifelse(length(llm_genes)==0, NA_real_,
                         length(intersect(llm_genes, gt_top20)) / length(llm_genes)),
      recall = length(intersect(llm_genes, gt_top20)) / length(gt_top20),
      jaccard_vs_truth = jaccard_set(llm_genes, gt_top20),
      overlap = overlap_coef(llm_genes, gt_top20),
      exact_match = setequal(llm_genes, gt_top20),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(scores)
}

pairwise_overlap_runs <- function(folder, llm_name) {
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 2)
  
  sets <- lapply(files, read_top20, key="top20_genes")
  names(sets) <- basename(files)
  
  n <- length(sets)
  O <- matrix(NA_real_, n, n, dimnames=list(names(sets), names(sets)))
  
  for (i in 1:n) {
    for (j in 1:n) {
      O[i,j] <- overlap_coef(sets[[i]], sets[[j]])
    }
  }
  
  vals <- O[upper.tri(O)]
  
  data.frame(
    llm = llm_name,
    mean_pairwise_overlap = mean(vals, na.rm=TRUE),
    sd_pairwise_overlap = sd(vals, na.rm=TRUE),
    min_pairwise_overlap = min(vals, na.rm=TRUE),
    max_pairwise_overlap = max(vals, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
}

####Ground truth for Prompt 6 (0.05 < FDR ≤ 0.15)

gt_sig_p6 <- sort(as.character(gt$gene[!is.na(gt$padj) & gt$padj > 0.05 & gt$padj <= 0.15]))
gt_n_p6   <- length(gt_sig_p6)
cat("Ground-truth genes for Prompt 6 (0.05 < FDR ≤ 0.15):", gt_n_p6, "\n")

###################################################
########Run Prompt P6 for all 3 LLMs################
###################################################

# --- Change these paths if needed ---
p6_chatgpt_dir <- "ChatGPT/prompt_v6"
p6_gemini_dir  <- "Gemini/prompt_v6"
p6_claude_dir  <- "Claude/prompt_v6"

scores_p6_chatgpt <- score_folder_top20(p6_chatgpt_dir, "ChatGPT", gt_top20_p6)
scores_p6_gemini  <- score_folder_top20(p6_gemini_dir,  "Gemini",  gt_top20_p6)
scores_p6_claude  <- score_folder_top20(p6_claude_dir,  "Claude",  gt_top20_p6)

scores_p6_all <- dplyr::bind_rows(scores_p6_chatgpt, scores_p6_gemini, scores_p6_claude)

print(scores_p6_all)

summary_p6 <- scores_p6_all %>%
  group_by(llm) %>%
  summarise(
    runs = n(),
    mean_jaccard_vs_truth = mean(jaccard_vs_truth, na.rm=TRUE),
    exact_match_rate = mean(exact_match, na.rm=TRUE),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall = mean(recall, na.rm=TRUE),
    mean_overlap = mean(overlap, na.rm=TRUE)
  )


print(summary_p6)

#dir.create("results", showWarnings = FALSE)
#write.csv(scores_p6_all, "results/prompt6_run_level_scores.csv", row.names=FALSE)
#write.csv(summary_p6,   "results/prompt6_summary_by_llm.csv", row.names=FALSE)



# Truth set for prompt v6: borderline window, then deterministic Top-20 rule
# Rule: smallest padj, break ties by largest |log2FoldChange|
gt_top20_p6 <- gt %>%
  filter(!is.na(padj), padj > 0.05, padj <= 0.15) %>%
  mutate(abs_lfc = abs(log2FoldChange)) %>%
  arrange(padj, desc(abs_lfc)) %>%
  slice_head(n = 20) %>%
  pull(gene) %>%
  as.character()

gt_top20_p6 <- sort(unique(gt_top20_p6))
gt_n_top20  <- length(gt_top20_p6)

cat("Truth Top-20 size (prompt v6):", gt_n_top20, "\n")
print(gt_top20_p6)


########Within llm stability 

stab_p6 <- dplyr::bind_rows(
  pairwise_jaccard_runs(p6_chatgpt_dir, "ChatGPT"),
  pairwise_jaccard_runs(p6_gemini_dir,  "Gemini"),
  pairwise_jaccard_runs(p6_claude_dir,  "Claude")
)

print(stab_p6)


stab_overlap_p6 <- dplyr::bind_rows(
  pairwise_overlap_runs(p6_chatgpt_dir, "ChatGPT"),
  pairwise_overlap_runs(p6_gemini_dir, "Gemini"),
  pairwise_overlap_runs(p6_claude_dir, "Claude")
)

print(stab_overlap_p6)

#dir.create("../results", showWarnings = FALSE)   # if you're inside results/llm_tables

#write.csv(scores_p6_all, "../results/prompt_v6_run_level_scores.csv", row.names=FALSE)
#write.csv(summary_p6,    "../results/prompt_v6_summary_by_llm.csv", row.names=FALSE)
#write.csv(stab_p6,       "../results/prompt_v6_within_llm_stability.csv", row.names=FALSE)

###################################
######Accuracy vs Stability Scatter
df_plot <- summary_p6 %>%
  left_join(stab_p6, by = "llm")

ggplot(df_plot, aes(
  x = mean_pairwise_jaccard,
  y = mean_jaccard_vs_truth,
  label = llm
)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8) +
  xlim(0.5,1.5) + ylim(0,1) +
  theme(
    axis.title = element_text(size=12),
    axis.text = element_text(size=11)
  ) +
  labs(
    title = "Stability does not imply correctness",
    x = "Within-LLM stability (Jaccard)",
    y = "Agreement with ground truth"
  )

df_plot_overlap <- summary_p6 %>%
  left_join(stab_overlap_p6, by = "llm")

ggplot(df_plot_overlap, aes(
  x = mean_pairwise_overlap,
  y = mean_jaccard_vs_truth,
  label = llm
)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8) +
  xlim(0.5,1.5) + ylim(0,1) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  labs(
    title = "Stability vs correctness (Overlap)",
    x = "Within-LLM stability (Overlap coefficient)",
    y = "Agreement with ground truth"
  )


##################################################
############D: promnpt P7a verus P7b############
#################################################

library(jsonlite)
library(dplyr)

jaccard_set <- function(a, b) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

overlap_coef <- function(a, b) {
  a <- unique(as.character(a)); b <- unique(as.character(b))
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / min(length(a), length(b))
}


read_top20 <- function(f, key="top20_genes") {
  out <- fromJSON(f)
  genes <- out[[key]]
  if (is.null(genes)) genes <- character(0)
  if (is.list(genes)) genes <- unlist(genes, use.names = FALSE)
  sort(unique(as.character(genes)))
}

# Pairwise Jaccard across runs inside one folder
within_prompt_stability <- function(folder, llm_name, prompt_label) {
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 2)
  
  sets <- lapply(files, read_top20)
  n <- length(sets)
  
  J <- matrix(NA_real_, n, n)
  O <- matrix(NA_real_, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      J[i,j] <- jaccard_set(sets[[i]], sets[[j]])
      O[i,j] <- overlap_coef(sets[[i]], sets[[j]])
    }
  }
  
  j_vals <- J[upper.tri(J)]
  o_vals <- O[upper.tri(O)]
  
  data.frame(
    llm = llm_name,
    prompt = prompt_label,
    runs = n,
    mean_pairwise_jaccard = mean(j_vals, na.rm=TRUE),
    sd_pairwise_jaccard = sd(j_vals, na.rm=TRUE),
    mean_pairwise_overlap = mean(o_vals, na.rm=TRUE),
    sd_pairwise_overlap = sd(o_vals, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
}

# Canonical set for a prompt:
# - strict: intersection across runs (very conservative)
canonical_intersection <- function(folder) {
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 1)
  Reduce(intersect, lapply(files, read_top20))
}

# - majority vote: include genes appearing in >= ceil(runs/2)
canonical_majority <- function(folder) {
  files <- list.files(folder, pattern="\\.json$", full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(files) >= 1)
  sets <- lapply(files, read_top20)
  allg <- unlist(sets, use.names = FALSE)
  tab <- table(allg)
  thr <- ceiling(length(sets) / 2)
  sort(names(tab)[tab >= thr])
}


chatgpt_7a <- "ChatGPT/prompt_v7a"
chatgpt_7b <- "ChatGPT/prompt_v7b"

gemini_7a  <- "Gemini/prompt_v7a"
gemini_7b  <- "Gemini/prompt_v7b"

claude_7a  <- "Claude/prompt_v7a"
claude_7b  <- "Claude/prompt_v7b"

#####Within-prompt stability (paper table)

stab_7 <- bind_rows(
  within_prompt_stability(chatgpt_7a, "ChatGPT", "7a"),
  within_prompt_stability(chatgpt_7b, "ChatGPT", "7b"),
  within_prompt_stability(gemini_7a,  "Gemini",  "7a"),
  within_prompt_stability(gemini_7b,  "Gemini",  "7b"),
  within_prompt_stability(claude_7a,  "Claude",  "7a"),
  within_prompt_stability(claude_7b,  "Claude",  "7b")
)

print(stab_7)

######################################################################
####////// Prompt sensitivity within each LLM (P7a vs P7b)\\\\\\######
######################################################################

canon_7 <- function(dir7a, dir7b, llm_name) {
  
  A <- canonical_majority(dir7a)
  B <- canonical_majority(dir7b)
  
  data.frame(
    llm = llm_name,
    jaccard_7a_vs_7b = jaccard_set(A, B),
    overlap_7a_vs_7b = overlap_coef(A, B),   # ← ADD THIS LINE
    n_7a = length(A),
    n_7b = length(B),
    n_overlap = length(intersect(A, B)),
    n_only_7a = length(setdiff(A, B)),
    n_only_7b = length(setdiff(B, A)),
    stringsAsFactors = FALSE
  )
}
sens_7 <- bind_rows(
  canon_7(chatgpt_7a, chatgpt_7b, "ChatGPT"),
  canon_7(gemini_7a,  gemini_7b,  "Gemini"),
  canon_7(claude_7a,  claude_7b,  "Claude")
)

print(sens_7)

export_flips <- function(dir7a, dir7b, llm_name) {
  A <- canonical_majority(dir7a)
  B <- canonical_majority(dir7b)
  
  data.frame(
    llm = llm_name,
    only_in_7a = paste(setdiff(A, B), collapse = ";"),
    only_in_7b = paste(setdiff(B, A), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

flips_7 <- bind_rows(
  export_flips(chatgpt_7a, chatgpt_7b, "ChatGPT"),
  export_flips(gemini_7a,  gemini_7b,  "Gemini"),
  export_flips(claude_7a,  claude_7b,  "Claude")
)

print(flips_7)

chatgpt_7a <- "ChatGPT/prompt_v7a"
chatgpt_7b <- "ChatGPT/prompt_v7b"

A <- canonical_majority(chatgpt_7a); A <- head(A, 20)
B <- canonical_majority(chatgpt_7b); B <- head(B, 20)


#Prompt Sensitivity Plot (7a vs 7b)
ggplot(sens_7, aes(x = llm, y = jaccard_7a_vs_7b, fill = llm)) +
  geom_col(width = 0.6) +
  ylim(0,1) +
  theme(
    axis.title = element_text(size=12),
    axis.text = element_text(size=11)
  ) +
  labs(
    title = "Prompt sensitivity: 7a vs 7b",
    x = "LLM",
    y = "Jaccard similarity (Top-20 sets)"
  ) +
  guides(fill = "none")

ggplot(sens_7, aes(x = llm, y = overlap_7a_vs_7b, fill = llm)) +
  geom_col(width = 0.6) +
  ylim(0,1) +
  labs(
    title = "Prompt sensitivity (Overlap coefficient)",
    x = "LLM",
    y = "Overlap coefficient"
  ) +
  theme_minimal() +
  guides(fill = "none")

library(cowplot)
library(tidyr)

sens_long <- sens_7 |>
  pivot_longer(
    cols = c(jaccard_7a_vs_7b, overlap_7a_vs_7b),
    names_to = "metric",
    values_to = "value"
  )

ggplot(sens_long, aes(x = llm, y = value, fill = metric)) +
  geom_col(position = "dodge", width = 0.6) +
  ylim(0,1) +
  labs(
    title = "Prompt sensitivity (7a vs 7b)",
    x = "LLM",
    y = "Similarity (Top-20 gene sets)",
    fill = "Metric"
  ) +
  scale_fill_manual(
    values = c("#E64B35", "#4DBBD5"),
    labels = c("Jaccard similarity", "Overlap coefficient")
  ) +
  theme_minimal()

#########################################
########////E: promnpt P9\\\\##########
#########################################

library(jsonlite)
library(dplyr)

# You said wd is results/llm_tables
gt <- read.csv("DESeq2_top200_minimal.csv", stringsAsFactors = FALSE)

valid_genes <- sort(unique(as.character(gt$gene)))
length(valid_genes)

######Helper
library(jsonlite)

load_ranked_top20 <- function(f) {
  out <- fromJSON(f)
  rt <- out$ranked_top20
  if (is.null(rt)) stop("Missing ranked_top20 in ", f)
  
  if (is.data.frame(rt)) {
    genes <- as.character(rt$gene)
    ranks <- as.integer(rt$rank)
  } else {
    genes <- as.character(sapply(rt, `[[`, "gene"))
    ranks <- as.integer(sapply(rt, `[[`, "rank"))
  }
  
  if (length(genes) != 20)
    stop("Not 20 genes in ", f)
  
  if (!all(sort(ranks) == 1:20))
    stop("Ranks not 1..20 in ", f)
  
  if (anyDuplicated(genes))
    stop("Duplicate genes in ", f)
  
  setNames(ranks, genes)
}


############################################
####//Validity check for hallucinations\\##
############################################

valid_genes <- sort(unique(as.character(gt$gene)))

check_validity <- function(folder, llm_name) {
  files <- list.files(folder, "\\.json$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    ranks <- load_ranked_top20(f)
    invalid <- setdiff(names(ranks), valid_genes)
    
    data.frame(
      llm = llm_name,
      file = basename(f),
      n_invalid = length(invalid),
      invalid_genes = ifelse(length(invalid) == 0, "",
                             paste(invalid, collapse = ";")),
      stringsAsFactors = FALSE
    )
  }))
}

p9_valid <- rbind(
  check_validity("ChatGPT/prompt_v9", "ChatGPT"),
  check_validity("Gemini/prompt_v9",  "Gemini"),
  check_validity("Claude/prompt_v9",  "Claude")
)
print(p9_valid)

library(ggplot2)

######bar plot of invalid by llms
ggplot(p9_valid, aes(x = llm, y = n_invalid, fill = llm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme(
    axis.title = element_text(size=12),
    axis.text = element_text(size=11)
  ) +
  labs(
    title = "Hallucinated genes by LLM (Prompt 9)",
    x = "LLM",
    y = "Number of invalid genes per run"
  ) +
  guides(fill = "none")

#Narrative payoff:
#Claude is internally stable but consistently hallucinates, whereas ChatGPT and Gemini remain grounded


#Within-LLM Stability Heatmaps (Rank Stability)
#Already computed Spearman / Kendall matrices (S, K) earlier.
#Plot: Spearman correlation heatmap (per LLM, per prommpt)


########Within llm rank stability 

spearman_cor <- function(x,y) suppressWarnings(cor(x,y,method="spearman"))
kendall_cor  <- function(x,y) suppressWarnings(cor(x,y,method="kendall"))
jaccard_set  <- function(a,b) length(intersect(a,b))/length(union(a,b))

rank_stability <- function(folder, llm_name) {
  files <- list.files(folder, "\\.json$", full.names = TRUE)
  if (length(files) < 2) stop("Need ≥2 runs in ", folder)
  maps <- lapply(files, load_ranked_top20)
  n <- length(maps)
  
  S <- K <- J <- matrix(NA, n, n)
  
  for (i in 1:n) for (j in 1:n) {
    common <- intersect(names(maps[[i]]), names(maps[[j]]))
    S[i,j] <- spearman_cor(maps[[i]][common], maps[[j]][common])
    K[i,j] <- kendall_cor(maps[[i]][common],  maps[[j]][common])
    J[i,j] <- jaccard_set(names(maps[[i]]), names(maps[[j]]))
  }
  data.frame(
    llm = llm_name,
    mean_spearman = mean(S[upper.tri(S)], na.rm=TRUE),
    mean_kendall  = mean(K[upper.tri(K)], na.rm=TRUE),
    mean_jaccard  = mean(J[upper.tri(J)], na.rm=TRUE)
  )
}

p9_stability <- rbind(
  rank_stability("ChatGPT/prompt_v9", "ChatGPT"),
  rank_stability("Gemini/prompt_v9",  "Gemini"),
  rank_stability("Claude/prompt_v9",  "Claude")
)
print(p9_stability)

####Across-Prompt Summary (Optional but Clean)
#Aggregate per LLM across selected prompts (1, 5, 6, 7, 9):
#  LLM	Prompt	Stability	Validity	Prompt Sensitivity

all_prompt_summary <- bind_rows(
  summary_p1 %>% mutate(prompt = "P1 (FDR≤0.05)"),
  summary_p5 %>% mutate(prompt = "P5 (FDR≤0.10)"),
  summary_p6 %>% mutate(prompt = "P6 (Borderline Top-20)"),
  sens_7     %>% mutate(prompt = "P7 (7a vs 7b)") %>%
    rename(mean_jaccard = jaccard_7a_vs_7b)
)

ggplot(all_prompt_summary,
       aes(x = prompt, y = mean_jaccard, color = llm, group = llm)) +
  geom_line() +
  geom_point(size = 3) +
  theme(
    axis.title = element_text(size=12),
    axis.text = element_text(size=11)
  ) +
  labs(
    title = "LLM behavior across prompt regimes",
    y = "Mean Jaccard similarity",
    x = "Prompt"
  )

ggplot(all_prompt_summary,
       aes(x = prompt, y = mean_overlap, color = llm, group = llm)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    title = "LLM behavior across prompt regimes (Overlap coefficient)",
    y = "Overlap coefficient",
    x = "Prompt"
  ) +
  theme_minimal()


