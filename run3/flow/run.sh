#Parameters
#----------

#- config (str): path of directory with config files
#- an_res_file (str): path of directory with analysis results
#- centrality (str): centrality class
#- resolution (str/int): resolution file or resolution value
#- outputdir (str): output directory
#- suffix (str): suffix for output files
#- vn_method (str): vn technique (sp, ep, deltaphi)
#- wagon_id (str): wagon ID
#- skip_resolution (bool): skip resolution extraction
#- skip_projection (bool): skip projection extraction
#- skip_vn (bool): skip raw yield extraction
#----------

# resutls will be saved in ${workdir}/Results/${dataCent}/${centrality}/${size}/
scriptdir=$(dirname $0)
workdir=/home/wuct/ALICE/local/DmesonAnalysis
config=/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config/config_flow_RESO.yml

centrality=k3050 # k020 k3050 k6080  <-----------------------------------------------------------------------------------------------
dataCent=2060 # 020 2050 50100   <---------------------------------------------------------------------------------------------------
size=large # small medium large <-------------------------------------------------------------------------------------------------------
vn_method=sp # sp ep deltaphi   <---------------------------------------------------------------------------------------------------
qvec= # full recenter   <-------------------------------------------------------------------------------------------------------
debug= #_old #_new #_tot

wagon_id= # 13649 14351 13650 14352    <---------------------------------------------------------------------------------------
doReso=true # false true    <-------------------------------------------------------------------------------------------------------
doProj=false # false true

an_res_file="/media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_0.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_1.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_2.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_3.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_4.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_5.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_6.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_7.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_8.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_9.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_10.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_11.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_12.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_13.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_14.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_15.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_16.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_17.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_18.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_19.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_20.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_21.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_22.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_23.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_24.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_25.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_26.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_27.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_28.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_29.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_30.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_31.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_32.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_33.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_34.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_35.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_36.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_37.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_38.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_39.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_40.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_41.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_42.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_43.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_44.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_45.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_46.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_47.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_48.root \ 
  /media/wuct/wulby/ALICE/AnRes/D0_flow/2024/Reso/389979/temp_merged_s0_49.root
"

resolution=path/to/resolution.root

if [ ! -z "$wagon_id" ]; then wagon="-w ${wagon_id}" ; else wagon="" ; fi

# suffix
cent="${centrality:1}"
siz="${size:0:1}"
qve="${qvec:0:2}"
suffix=${cent}${siz}_${qve}${debug}

if $doReso; then
    reso="--skip_projection --skip_vn"
    suffix=${suffix}_Reso
else 
    reso="-r ${resolution} --skip_resolution" 
fi

if $doProj; then
    proj="--skip_resolution --skip_vn"
    suffix=${suffix}_Proj
else
    proj=""
fi

# output dir.
outputdir=${workdir}/Results/${dataCent}/${centrality}/${size}/
if [ ! -d "${outputdir}" ]; then mkdir -p ${outputdir}; fi

    python3 ${scriptdir}/run_full_flow_analysis.py \
    ${config} \
    ${an_res_file} \
    -c ${centrality} \
    -o ${outputdir} \
    -s ${suffix} \
    -v ${vn_method} \
    ${reso} \
    ${wagon} \
    $proj \
    --skip_efficiency \
    #--skip_projection

if [ ! -z "$wagon_id" ]; then
    cp -rf ${outputdir}/${wagon_id}/*  ${outputdir}
    rm -rf ${outputdir}/${wagon_id}
fi

