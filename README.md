# LLM Gene Selection Prompt Sensitivity

This repository contains code and evaluation outputs for the paper: "Stability Does Not Imply Correctness: Prompt Sensitivity in LLM-based Gene Selection."


This repository contains prompts, raw LLM outputs, and evaluation scripts used to analyze prompt sensitivity in LLM-based gene prioritization tasks using a fixed DESeq2 reference.

## Repository Structure

data/
DESeq2 reference tables used as statistical ground truth.

prompts/
All prompt templates used in the experiments.

outputs/
Raw LLM responses for each model and prompt configuration.

scripts/
R scripts used to compute Jaccard similarity and overlap coefficients.

figures/
Plots included in the paper.

## Models evaluated

ChatGPT (GPT-5.2)  
Gemini 3  
Claude Opus 4.5

## Reproducibility

## Quick start

1. Install R (>=4.2)
2. Install required packages:

install.packages(c("jsonlite","tidyverse","data.table"))

3. Run:

Rscript scripts/run_analysis.R

Running `scripts/run_analysis.R` will recompute all metrics and figures from the raw outputs.

## License

MIT License
