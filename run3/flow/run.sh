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

workdir=/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/
config=/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/Results/2060/k3050/large/sp/config/config_flow_3050l_PASS4_full_PbPb.yml
#config=/home/wuct/localAnalysis/flow/DmesonAnalysis/run3/flow/Results/2060/k3050/large/sp/config/config_flow_3050l__set3_pol2_2_16.yml

centrality=k3050 # k020 k3050 k6080  <-----------------------------------------------------------------------------------------------
dataCent=2060 # 020 2060 50100   <---------------------------------------------------------------------------------------------------
size=large # small medium large <-------------------------------------------------------------------------------------------------------
vn_method=sp # sp ep deltaphi   <---------------------------------------------------------------------------------------------------
qvec= # full recenter   <-------------------------------------------------------------------------------------------------------
debug=PASS4_full_PbPb_0 #_old #_rossi_3050_klin #_254365_2_24 #_new #_tot_old_allbins #_tot _Largebin

wagon_id= # 13649 14351 13650 14352    <---------------------------------------------------------------------------------------
doReso=false # false true    <-------------------------------------------------------------------------------------------------------
doProj=false # false true

# an_res_file="/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/k6080/full/temp_merged_s2_0.root \
# /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/k6080/full/temp_merged_s2_1.root \
# /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/k6080/full/temp_merged_s2_2.root \
# /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/k6080/full/temp_merged_s2_3.root"
an_res_file="/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_0.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_1.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_2.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_3.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_4.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_5.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_6.root \
/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s3_7.root
"
# an_res_file=/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/CombPID_3SigCut_AND/AnalysisResults_data_goldenRun.root
resolution=/media/wuct/wulby/ALICE/AnRes/resolution/output_reso/resosp3050l_PASS4_full_PbPb_Reso.root

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

    python3 ${workdir}/run_full_flow_analysis.py \
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

