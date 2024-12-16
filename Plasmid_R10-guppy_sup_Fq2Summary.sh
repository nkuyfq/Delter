#!/usr/bin/env bash

Rscript scripts/Plasmid_R10_TP_FP_Qscore.R ${snakemake_params[opts1]} "${snakemake_params[opts2]}" "${snakemake_output[0]}"
