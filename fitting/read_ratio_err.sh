#! /bin/bash

TODAY="20140223"
SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError"

DIRECTORY="${TODAY}_${SUFFIX}"


bkgfunctions=(pol1 pol2 pol3 pol4 pol5 pol6 expFunct gausBkg)
ptbins=(65-30 3-30 3-65)
rapbins=(0-24 0-16 16-24)
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
	if [[ ${rap} != "16-24" && ${pt} != "65-30" ]];
	then
	    continue;
	fi;
	for cent in "${centbins[@]}";
	do
	    if [[ ${cent} != "0-100" && ${pt} != "3-30" && ${rap} == "16-24" ]];
	    then
		continue;
	    fi;
	    for bkg in "${bkgfunctions[@]}";
	    do
		echo ${bkg} ${ratio_tmp} >> fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
		if [[ ((${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol1" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol1" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "16-24" && ${bkg} == "pol3" && ${cent} == "0-100") ||
		       (${pt} == "3-30" && ${rap} == "16-24" && ${bkg} == "pol3" && ${cent} == "0-100") ||
		       (${pt} == "3-65"  && ${rap} == "16-24" && ${bkg} == "pol2" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol1" && ${cent} == "0-20")  ||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol1" && ${cent} == "0-20")  ||
		       (${pt} == "3-30" && ${rap} == "16-24" && ${bkg} == "pol3" && ${cent} == "0-20")  ||
		       (${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol2" && ${cent} == "20-40") ||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol1" && ${cent} == "20-40") ||
		       (${pt} == "3-30" && ${rap} == "16-24" && ${bkg} == "pol2" && ${cent} == "20-40") ||
		       (${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol1" && ${cent} == "40-100")||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol1" && ${cent} == "40-100")||
		       (${pt} == "3-30" && ${rap} == "16-24" && ${bkg} == "pol1" && ${cent} == "40-100")
		      ) ]];
		then
		ratio=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
		ratio_err=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $5}'`
		echo "PbPb rap" ${rap} " pT"${pt} " cent" ${cent} ${bkg}": "
		echo ${ratio}
		echo ${ratio_err}
		echo "---"
		fi;
	    done;
	done;
    done;
done;


for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${rap} != "16-24" && ${pt} != "65-30" ]];
	then
	    continue;
	fi;
	for cent in "0-100";
	do
	    for bkg in "${bkgfunctions[@]}";
	    do
		ratio_tmp=`grep 'fracP =' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
		echo ${bkg} ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
		if [[ ((${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol3" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "16-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3-30"  && ${rap} == "16-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3-65"  && ${rap} == "16-24" && ${bkg} == "pol1" && ${cent} == "0-100")
		      ) ]];
		then
		    ratio=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
		    ratio_err=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $5}'`
		    echo "pp rap" ${rap} " pT"${pt} " cent" ${cent} ${bkg}": "
		    echo ${ratio}
		    echo ${ratio_err}
		    echo "---"
		fi;
	    done;
	done;
    done;
done;
