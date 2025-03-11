#!/bin/bash
export V2VsFracFile="/home/wuct/ALICE/local/Results/BDT/k6080/third/syst/pre_sys/cutvar_combined/V2VsFrac/V2VsFrac_combined.root"
export config_syst="/home/wuct/ALICE/local/sys_dev/DmesonAnalysis/run3/flow/systematics/syst_summary.yml"
export Dmeson="D0"
export cent="k6080"
export suffix="syst"
export output_dir="/home/wuct/ALICE/local/Results/BDT/k6080/third/syst/syst_total"

python3 compute_promptvn_withsyst_ste.py $V2VsFracFile $config_syst -d $Dmeson -c $cent -s $suffix -o $output_dir

