import yaml
import os
mc_files = [
    "/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/flow/MC/AnalysisResults_full_default_407164.root",
    "/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/flow/MC/AnalysisResults_full_ptsmearing1p5_407162.root",
    "/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/flow/MC/AnalysisResults_full_ptsmearing1p5_vsphi_407161.root",
]
track_tuner_lables = [
    "default",
    "ptsmearing1p5",
    "ptsmearing1p5_vsphi"
]

# track_tuner_lables = [
#     "default",
#     "ptsmearing",
#     "ptsmearing_vsphi",
#     "ptsmearing1p5",
#     "ptsmearing1p5_vsphi",
#     "vsphi"
# ]

addtional_lables = "full"

config_file = "/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config/2024/k3050/config_flow_2024_small.yml"

for iFile, mc_file in enumerate(mc_files):
    print(f'Processing {mc_file}')
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    config['eff_filename'] = mc_file
    config['suffix'] = f"{track_tuner_lables[iFile]}_{addtional_lables}"
    with open(config_file, 'w') as file:
        yaml.dump(config, file)
    print('Config file updated:', config_file)
    os.system("cd /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT |bash /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.sh")

eff_files = []
for iFile, mc_file in enumerate(mc_files):
    eff_files.append(f'/home/wuct/ALICE/local/Results/BDT/2024/k3050/eff/cutvar_{track_tuner_lables[iFile]}_{addtional_lables}/eff/eff_{track_tuner_lables[iFile]}_{addtional_lables}_00.root')
    eff_files.append(f'/home/wuct/ALICE/local/Results/BDT/2024/k3050/eff/cutvar_{track_tuner_lables[iFile]}_{addtional_lables}/eff/eff_{track_tuner_lables[iFile]}_{addtional_lables}_00.root')
    eff_files
config_file_compare = "/home/wuct/ALICE/local/reso/DmesonAnalysis/comparisons/config_comparison_eff_tracktuner.yml"
with open(config_file_compare, 'r') as file:
    config = yaml.safe_load(file)
    config['inputs']['filenames'] = eff_files
    config['output']['filename'] = f"/home/wuct/ALICE/local/Results/BDT/2024/k3050/eff/pass4_eff_{addtional_lables}.root"
    with open(config_file_compare, 'w') as file:
        yaml.dump(config, file)
    print('Config file updated:', config_file_compare)
# os.system("cd /home/wuct/ALICE/local/reso/DmesonAnalysis/comparisons | python3 /home/wuct/ALICE/local/reso/DmesonAnalysis/comparisons/CompareGraphs.py /home/wuct/ALICE/local/reso/DmesonAnalysis/comparisons/config_comparison_eff_tracktuner.yml")

    