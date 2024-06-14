#!/usr/bin/env bash

export HDF5_PLUGIN_PATH=lib/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin

#bash scripts/delinfo2Signal+Qinfo_workflow.sh "${snakemake_input[0]}" "${snakemake_input[1]}" "${snakemake_input[2]}" "${snakemake_input[3]}" "${snakemake_input[4]}" "${snakemake_input[5]}" "${snakemake_input[6]}" "${snakemake_params[opts1]}" "${snakemake_params[opts2]}" "${snakemake_params[opts3]}" "${snakemake_params[opts4]}" "${snakemake_params[opts5]}" "${snakemake_output[0]}" "${snakemake_output[1]}" "${snakemake_output[2]}"

#bash scripts/delinfo2Signal+Qinfo_workflow-v2.sh "${snakemake_input[0]}" "${snakemake_input[1]}" "${snakemake_input[2]}" "${snakemake_input[3]}" "${snakemake_input[4]}" "${snakemake_input[5]}" "${snakemake_input[6]}" "${snakemake_params[opts1]}" "${snakemake_params[opts2]}" "${snakemake_params[opts3]}" "${snakemake_params[opts4]}" "${snakemake_params[opts5]}" "${snakemake_output[0]}" "${snakemake_output[1]}" "${snakemake_output[2]}"
bash scripts/delinfo2Signal+Qinfo_workflow-v3.sh "${snakemake_input[0]}" "${snakemake_input[1]}" "${snakemake_input[2]}" "${snakemake_input[3]}" "${snakemake_input[4]}" "${snakemake_input[5]}" "${snakemake_input[6]}" "${snakemake_params[opts1]}" "${snakemake_params[opts2]}" "${snakemake_params[opts3]}" "${snakemake_params[opts4]}" "${snakemake_params[opts5]}" "${snakemake_output[0]}" "${snakemake_output[1]}" "${snakemake_output[2]}"
