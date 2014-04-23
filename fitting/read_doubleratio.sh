#! /bin/bash

TODAY="20140411"
SUFFIX="SimFits_M2242_DblMu0_AllCent_2CB"

DIRECTORY="${TODAY}_${SUFFIX}"


#bkgfunctions=(pol0 pol1 pol2 pol3 pol4 pol5 pol6 expPol1 expPol2 expPol3 expPol4 expPol5 expPol6 expPol7)
ptbins=(65-30 3-30)
rapbins=(0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)

ratio_tmp=0;
ratio=0;
relRatio=0;
min=0;
max=0;

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ (${rap} == "16-24" && ${pt} == "65-30")|| 
	      (${rap} != "16-24" && ${pt} == "3-30")]];
	then
	    continue;
	fi;
	if [[ ${rap} == "0-24" ]];
	then
	    rapAlt="0.0-2.4"
	elif [[ ${rap} == "0-16" ]];
	then
	    rapAlt="0.0-1.6"
	elif [[ ${rap} == "16-24" ]];
	then
	    rapAlt="1.6-2.4"
	fi;	    
	if [[ ${pt} == "65-30" ]]
	then
	    ptAlt="6.5-30.0"
	elif [[ ${pt} == "3-30" ]]
	then
	    ptAlt="3.0-30.0"
	fi;
	echo ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt
	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI020 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI2040 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI40100 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI0100 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	zgrep -A4 'HI0100:' ${DIRECTORY}/sim_fitCBG_polBkg_rap${rapAlt}_pt${ptAlt}_cent.log.gz|grep "Double Ratio" |awk '{print $3; print $5}'
    done;
done;

echo "----------------------------"
echo "Signal shape"
echo "----------------------------"

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ (${rap} == "16-24" && ${pt} == "65-30")|| 
	      (${rap} != "16-24" && ${pt} == "3-30")]];
	then
	    continue;
	fi;
	if [[ ${rap} == "0-24" ]];
	then
	    rapAlt="0.0-2.4"
	elif [[ ${rap} == "0-16" ]];
	then
	    rapAlt="0.0-1.6"
	elif [[ ${rap} == "16-24" ]];
	then
	    rapAlt="1.6-2.4"
	fi;	    
	if [[ ${pt} == "65-30" ]]
	then
	    ptAlt="6.5-30.0"
	elif [[ ${pt} == "3-30" ]]
	then
	    ptAlt="3.0-30.0"
	fi;
	echo ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt
	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI020 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI2040 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI40100 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI0100 = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_centMult_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
	zgrep -A4 'HI0100:' ${DIRECTORY}/sim_fitCBG_polBkg_rap${rapAlt}_pt${ptAlt}_cent_allFree.log.gz|grep "Double Ratio" |awk '{print $3; print $5}'
    done;
done;


echo "----------------------------"
echo "Background shape"
echo "----------------------------"


for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ (${rap} == "16-24" && ${pt} == "65-30")|| 
	      (${rap} != "16-24" && ${pt} == "3-30")]];
	then
	    continue;
	fi;
	if [[ ${rap} == "0-24" ]];
	then
	    rapAlt="0.0-2.4"
	elif [[ ${rap} == "0-16" ]];
	then
	    rapAlt="0.0-1.6"
	elif [[ ${rap} == "16-24" ]];
	then
	    rapAlt="1.6-2.4"
	fi;	    
	if [[ ${pt} == "65-30" ]]
	then
	    ptAlt="6.5-30.0"
	elif [[ ${pt} == "3-30" ]]
	then
	    ptAlt="3.0-30.0"
	fi;
	echo ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt
	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI020 = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI2040 = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	grep 'doubleRatio_HI40100 = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI0100 = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_centMult_allVars.txt | awk '{print $3; print $5}'
	zgrep -A4 'HI0100:' ${DIRECTORY}/sim_fitCBG_expBkg_rap${rapAlt}_pt${ptAlt}_cent.log.gz|grep "Double Ratio" |awk '{print $3; print $5}'
    done;
done;


#TODAY="20140325"
#SUFFIX="SimFits_M2242_DblMu0_WithSyst_final"
#
#DIRECTORY="${TODAY}_${SUFFIX}"
#ptbins=(3-65 65-30)
#
#rap="16-24"
#for pt in "${ptbins[@]}";
#do
#	echo ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_allVars.txt
#	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_allVars.txt | awk '{print $3; print $5}'
#done;
#
#echo "----------------------------"
#echo "Signal shape"
#echo "----------------------------"
#
#for pt in "${ptbins[@]}";
#do
#	echo ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_freeAlpha_freeN_freeGwidth_allVars.txt
#	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI = ' ${DIRECTORY}/fracLogCBG_PbPbpol*_rap${rap}_pT${pt}_cent0-100_freeAlpha_freeN_freeGwidth_allVars.txt | awk '{print $3; print $5}'
#done;
#
#echo "----------------------------"
#echo "Background shape"
#echo "----------------------------"
#
#for pt in "${ptbins[@]}";
#do
#	echo ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_cent0-100_allVars.txt
#	grep 'fracP_pp = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_cent0-100_allVars.txt | awk '{print $3; print $5}'
#	grep 'doubleRatio_HI = ' ${DIRECTORY}/fracLogCBG_PbPbexpPol*_rap${rap}_pT${pt}_cent0-100_allVars.txt | awk '{print $3; print $5}'
#done;
