import os
import sys
import numpy as np
import argparse
import yaml
import shutil
sys.path.append('..')
from flow_analysis_utils import get_cut_sets_config
from ComputeDataDriFrac_flow import main_data_driven_frac
from ComputeV2vsFDFrac import main_v2_vs_frac
from concurrent.futures import ProcessPoolExecutor

def check_dir(dir):

    if not os.path.exists(dir):
        print(f"\033[32m{dir} does not exist, it will be created\033[0m")
        os.makedirs(dir)
    else:
        print(f"\033[33m{dir} already exists, it will be removed and recreat\033[0m")
        shutil.rmtree(dir)
        os.makedirs(dir)

    return

def process_proj_mc(i, ProjMcPath, config_flow, output_dir, suffix, ptweights_exists):
    iCutSets = f"{i:02d}"
    # TODO: load the path and bool from the config file
    given_weights = False
    if given_weights:
        ptweights_exists = True
        ptweightsPath = '/home/wuct/ALICE/local/Results/BDT/k3050/full/uncorrelated/cutvar_pt1_4/ptweights/pTweight_pt1_4.root'
    else:
        ptweightsPath = f'{output_dir}/ptweights/pTweight_{suffix}.root'
    if not ptweights_exists:
        print(f"\033[32mpython3 {ProjMcPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml -o {output_dir} -s {suffix}_{iCutSets}\033[0m")
        os.system(f"python3 {ProjMcPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml -o {output_dir} -s {suffix}_{iCutSets}")
    else:
        print(
            f"\033[32mpython3 {ProjMcPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml "
            f"-w {ptweightsPath} hPtWeightsFONLLtimesTAMUDcent "
            f"-wb {ptweightsPath} hPtWeightsFONLLtimesTAMUBcent "
            f"-o {output_dir} -s {suffix}_{iCutSets} \033[0m"
        )
        os.system(f"python3 {ProjMcPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml "
                  f"-w {ptweightsPath} hPtWeightsFONLLtimesTAMUDcent "
                  f"-wb {ptweightsPath} hPtWeightsFONLLtimesTAMUBcent "
                  f"-o {output_dir} -s {suffix}_{iCutSets}")

def process_efficiency_cutset(i, EffPath, config_flow, output_dir, suffix, cent):
    iCutSets = f"{i:02d}"
    print(f"\033[32mpython3 {EffPath} {config_flow} {output_dir}/proj_mc/proj_mc_{suffix}_{iCutSets}.root -c {cent} -o {output_dir} -s {suffix}_{iCutSets}\033[0m")
    print(f"\033[32mProcessing cutset {iCutSets}\033[0m")
    os.system(f"python3 {EffPath} {config_flow} {output_dir}/proj_mc/proj_mc_{suffix}_{iCutSets}.root -c {cent} -o {output_dir} -s {suffix}_{iCutSets} --batch")

def process_vn_cutset(i, SimFitPath, config_flow, cent, output_dir, suffix, vn_method):
    iCutSets = f"{i:02d}"
    print(f"\033[32mpython3 {SimFitPath} {config_flow} {cent} {output_dir}/proj/proj_{suffix}.root -o {output_dir}/ry -s _{suffix}_{iCutSets} -vn {vn_method}\033[0m")
    print(f"\033[32mProcessing cutset {iCutSets}\033[0m")
    os.system(f"python3 {SimFitPath} {config_flow} {cent} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -o {output_dir}/ry -s _{suffix}_{iCutSets} -vn {vn_method} --batch")

def run_full_cut_variation(config_flow, anres_dir, cent, res_file, output, suffix, vn_method, 
                           skip_calc_weights=False,
                           skip_make_yaml=False, 
                           skip_cut_variation=False,
                           skip_proj_mc=False,
                           skip_efficiency=False,
                           skip_vn = False,
                           skip_frac_cut_var=False,
                           skip_data_driven_frac=False,
                           skip_v2_vs_frac=False):

#___________________________________________________________________________________________________________________________
    # Load and copy the configuration file
    with open(config_flow, 'r') as cfgFlow:
        config = yaml.safe_load(cfgFlow)
    
    CutSets, _, _, _, _ = get_cut_sets_config(config_flow)
    nCutSets = max(CutSets)

    print(f"\033[32mNumber of cutsets: {nCutSets}\033[0m")

    output_dir = f"{output}/cutvar_{suffix}"
 
    os.system(f"mkdir -p {output_dir}")

    # the pT weights histograms
    PtWeightsDHistoName = 'hPtWeightsFONLLtimesTAMUDcent'
    PtWeightsBHistoName = 'hPtWeightsFONLLtimesTAMUBcent'
 
    # copy the configuration file
    config_suffix = 1
    os.makedirs(f'{output_dir}/config_flow', exist_ok=True)
    while os.path.exists(f'{output_dir}/config_flow/config_flow_{suffix}_{config_suffix}.yml'):
        config_suffix = config_suffix + 1
    os.system(f'cp {config_flow} {output_dir}/config_flow/config_flow_{suffix}_{config_suffix}.yml')
 
    # backup the results into history
    file_to_check = f"{output_dir}/V2VsFrac/V2VsFrac_{suffix}.root"
    if os.path.exists(file_to_check):
        for sub_path in ['ry', 'CutVarFrac', 'V2VsFrac']:
            os.system(f"mkdir -p {output_dir}/history/{config_suffix}/{sub_path}")
            os.system(f"cp {output_dir}/{sub_path}/* {output_dir}/history/{config_suffix}/{sub_path}")
        os.system(f"cp {output_dir}/config_flow/config_flow_{suffix}_{config_suffix-1}.yml {output_dir}/history/{config_suffix}")

#___________________________________________________________________________________________________________________________
    # calculate the pT weights
    if not skip_calc_weights:
        check_dir(f"{output_dir}/ptweights")
        CalcWeiPath = "./ComputePtWeights.py"

        print(f"\033[32mpython3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}\033[0m")
        os.system(f"python3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}")
    else:
        print("\033[33mWARNING: Calculation of weights will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # make yaml file
    if not skip_make_yaml:
        check_dir(f"{output_dir}/config")
        MakeyamlPath = './make_yaml_for_ml.py'

        print(f"\033[32mpython3 {MakeyamlPath} {config_flow} -o {output_dir} -s {suffix}\033[0m")
        os.system(f"python3 {MakeyamlPath} {config_flow} -o {output_dir} -s {suffix}")
    else:
        print("\033[33mWARNING: Make yaml will not be performed\033[0m")
    #TODO: 1.keep the yaml file for the user to check 2.modify the proj_thn_mc 3.use make_combination in proj_thn_mc.py

#___________________________________________________________________________________________________________________________
    # Cut variation (aply the cut and project)
    if not skip_cut_variation:
        check_dir(f"{output_dir}/proj")
        CutVarPath = "./cut_variation.py"

        anres_files = " ".join(anres_dir)
        print(f"\033[32mpython3 {CutVarPath} {config_flow} {anres_files} -c {cent} -r {res_file} -o {output_dir} -s {suffix}\033[0m")
        os.system(f"python3 {CutVarPath} {config_flow} {anres_files} -c {cent} -r {res_file} -o {output_dir} -s {suffix}")
    else:
        print("\033[33mWARNING: Cut variation will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # Projection for MC and apply the ptweights
    if not skip_proj_mc:
        check_dir(f"{output_dir}/proj_mc")
        ProjMcPath = "./proj_thn_mc.py"
        
        ptweights_exists = os.path.exists(f'{output_dir}/ptweights/pTweight_{suffix}.root')

        max_workers = 16 # hyper parameter default: 1
        with ProcessPoolExecutor(max_workers) as executor:
            futures = []
            for i in range(nCutSets):
                futures.append(executor.submit(process_proj_mc, i, ProjMcPath, config_flow, output_dir, suffix, ptweights_exists))
                sys.stdout.flush()
            
            for future in futures:
                future.result()

    else:
        print("\033[33mWARNING: Projection for MC will not be performed\033[0m")							

#___________________________________________________________________________________________________________________________
    # Compute the efficiency
    if not skip_efficiency:
        check_dir(f"{output_dir}/eff")
        EffPath = "./../compute_efficiency.py"

        max_workers = 32 # hyper parameter default: 1
        with ProcessPoolExecutor(max_workers) as executor:
            futures = []
            for i in range(nCutSets):
                futures.append(executor.submit(process_efficiency_cutset, i, EffPath, config_flow, output_dir, suffix, cent))
                sys.stdout.flush()
            
            for future in futures:
                future.result()
    else:
        print("\033[33mWARNING: Efficiency will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # do the simulation fit to get the raw yields
    if not skip_vn:
        check_dir(f"{output_dir}/ry")
        SimFitPath = "./../get_vn_vs_mass.py"
  
        max_workers = 32 # hyper parameter default: 1
        with ProcessPoolExecutor(max_workers) as executor:
            futures = []
            for i in range(nCutSets):
                futures.append(executor.submit(process_vn_cutset, i, SimFitPath, config_flow, cent, output_dir, suffix, vn_method))
                sys.stdout.flush()
            
            for future in futures:
                future.result()
    else:
        print("\033[33mWARNING: vn extraction will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # Compute the fraction by cut variation method
    if not skip_frac_cut_var:
        check_dir(f"{output_dir}/CutVarFrac")
        CurVarFracPath = "./compute_frac_cut_var.py"

        print(f"\033[32mpython3 {CurVarFracPath} {config_flow} {output_dir} -o {output_dir} -s {suffix}\033[0m")
        os.system(f"python3 {CurVarFracPath} {config_flow} {output_dir} -o {output_dir} -s {suffix} --batch")
    else:
        print("\033[33mWARNING: Fraction by cut variation will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # Compute fraction by Data-driven method
    if not skip_data_driven_frac:
        check_dir(f"{output_dir}/DataDrivenFrac")
        DataDrivenFracPath = "./ComputeDataDriFrac_flow.py"
        
        combined = config['combined'] if 'combined' in config else False
        print(f"\033[32mCombined method: {combined}\033[0m")
        if combined:
            if config['minimisation']['correlated']:
                print("\033[31mOnly support combined method when the minimisation is uncorrelated\033[0m")
                print("\033[31mThe combined method will not be performed\033[0m")
                # run the data-driven method with the corelated results
                print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b\033[0m")
                main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, combined=False)
            else:
                # the path of corresponding results with correlated cut method
                correlatedPath = config['correlatedPath'] if 'correlatedPath' in config else ''
                if correlatedPath == '':
                    print("\033[31mPlease provide the path to the corresponding correlated cut method\033[0m")
                    correlatedPath = input("Path(`output_dir` in run_cutvar.sh): ")
                correlatedCutVarPath = correlatedPath + '/cutvar_' + suffix
                
                # the path of the combined results
                outputdir_corr = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined')
                check_dir(f"{outputdir_corr}/DataDrivenFrac")
                
                # run the data-driven method with the combined results and the uncorrelated results
                print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b -comb -cc {correlatedCutVarPath} -oc {outputdir_corr}\033[0m") 
                main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, \
                                        combined=True, correlatedCutVarPath=correlatedCutVarPath, outputdir_combined=outputdir_corr)

        # run the data-driven method with the uncorrelated or correlated results
        else:
            print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b\033[0m")
            main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, combined=False)

    else:
        print("\033[33mWARNING: Fraction by Data-driven method will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
    # Compute v2 vs fraction
    if not skip_v2_vs_frac:
        check_dir(f"{output_dir}/V2VsFrac")
        v2vsFDFracPath = "./ComputeV2vsFDFrac.py"
        
        combined = config['combined'] if 'combined' in config else False
        
        if combined:
            inputdir_combined = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined')
            outputdir_combined = inputdir_combined
            check_dir(f"{outputdir_combined}/V2VsFrac")
            
            print(f"\033[32mthe combined method will be performed\033[0m")
            print(f"\033[32mpython3 {v2vsFDFracPath} {config_flow} -i {output_dir} -o {output_dir} -s {suffix} -comb -ic {inputdir_combined} -oc {outputdir_combined}\033[0m")
            main_v2_vs_frac(config=config_flow, inputdir=output_dir, outputdir=output_dir, suffix=suffix, \
                                combined=True, inputdir_combined=inputdir_combined, outputdir_combined=outputdir_combined)
        else:
            print(f"\033[32mpython3 {v2vsFDFracPath} {config_flow} -i {output_dir} -o {output_dir} -s {suffix}\033[0m")
            main_v2_vs_frac(config=config_flow, inputdir=output_dir, outputdir=output_dir, suffix=suffix, combined=False)
    else:
        print("\033[33mWARNING: v2 vs fraction will not be performed\033[0m")
    
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('flow_config', metavar='text', default='config_flow_d0.yml', help='configuration file')
    parser.add_argument('anres_dir', metavar='text', nargs='+', help='input ROOT files with anres')
    parser.add_argument("--centrality", "-c", metavar="text",default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r",  default="", help="resolution file/value")
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text", default="", help="suffix for output files")
    parser.add_argument("--vn_method", "-vn", metavar="text", default="sp", help="vn technique (sp, ep, deltaphi)")
    parser.add_argument("--skip_calc_weights", "-scw", action="store_true", help="skip calculation of weights")
    parser.add_argument("--skip_make_yaml", "-smy", action="store_true", help="skip make yaml")
    parser.add_argument("--skip_cut_variation", "-scv", action="store_true", help="skip cut variation")
    parser.add_argument("--skip_proj_mc", "-spm", action="store_true", help="skip projection for MC")
    parser.add_argument("--skip_efficiency", "-se", action="store_true", help="skip efficiency")
    parser.add_argument("--skip_vn", "-svn", action="store_true", help="skip vn extraction")
    parser.add_argument("--skip_frac_cut_var", "-sf", action="store_true", help="skip fraction by cut variation")
    parser.add_argument("--skip_data_driven_frac", "-sddf", action="store_true", help="skip fraction by data-driven method")
    parser.add_argument("--skip_v2_vs_frac", "-sv2fd", action="store_true", help="skip v2 vs FD fraction")
    args = parser.parse_args()

    run_full_cut_variation(args.flow_config, args.anres_dir, args.centrality, args.resolution, args.outputdir, args.suffix, args.vn_method, 
                        args.skip_calc_weights,
                        args.skip_make_yaml, 
                        args.skip_cut_variation, 
                        args.skip_proj_mc, 
                        args.skip_efficiency, 
                        args.skip_vn,
                        args.skip_frac_cut_var, 
                        args.skip_data_driven_frac, 
                        args.skip_v2_vs_frac)