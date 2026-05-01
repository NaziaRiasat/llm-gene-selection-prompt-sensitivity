# LLM Gene Selection Prompt Sensitivity

This repository contains code and evaluation outputs for the paper: "WHEN STABILITY FAILS: HIDDEN FAILURE MODES
OF LLMS IN DATA-CONSTRAINED SCIENTIFIC DECISION-MAKING."
Url: https://arxiv.org/pdf/2603.15840

This repository contains prompts, raw LLM outputs, and evaluation scripts used to analyze prompt sensitivity in LLM-based gene prioritization tasks using a fixed DESeq2 reference.

## Models evaluated

ChatGPT (GPT-5.2)  
Gemini 3  
Claude Opus 4.5


## Quick start

1. Install R (>=4.2)
2. Install required packages:

install.packages(c("jsonlite","tidyverse","data.table"))

3. Run:

Rscript scripts/run_analysis.R

Running `scripts/run_analysis.R` will recompute all metrics and figures from the raw outputs.

## Reproducing Results

1. Ground truth DESeq2 tables are in `data/`.
2. Prompt templates are in `prompts/`.
3. Raw LLM outputs are stored in `outputs/`.
4. Run the evaluation script:

   Rscript scripts/llm_score_output.R

This script computes:
- precision
- recall
- Jaccard similarity
- overlap coefficient
- exact match rate

and generates the summary statistics used in the paper.


## Data

Differential expression reference tables were derived from the NSCLC tumor-draining lymph node dataset: GEO accession: GSE239514

## Citation

If you use this code, please cite:

Code, prompts, raw LLM outputs, and evaluation scripts are available at:
https://github.com/NaziaRiasat/llm-prompt-sensitivity

WHEN STABILITY FAILS: HIDDEN FAILURE MODES OF LLMS IN DATA-CONSTRAINED SCIENTIFIC DECISION-MAKING
ICLR 2026 Workshop: I Can't Believe It's Not Better.

## License

MIT License
