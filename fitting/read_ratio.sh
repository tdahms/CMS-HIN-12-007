#! /bin/bash

TODAY="20140318"
SUFFIX="M2242_DblMu0"

DIRECTORY="${TODAY}_${SUFFIX}"


bkgfunctions=(pol0 pol1 pol2 pol3 pol4 pol5 pol6 expPol1 expPol2 expPol3 expPol4)
ptbins=(65-30 3-30 3-65)
rapbins=(0-24 0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)

ratio_tmp=0;
ratio=0;
relRatio=0;
min=0;
max=0;

for cent in "${centbins[@]}";
do
    for pt in "${ptbins[@]}";
    do
	for rap in "${rapbins[@]}";
	do
	    if [[ ${rap} != "16-24" && ${pt} != "65-30" ]];
	    then
		continue;
	    fi;
	    if [[ ${cent} != "0-100" && ${pt} != "3-30" && ${rap} == "16-24" ]];
	    then
		continue;
	    fi;
	    for bkg in "${bkgfunctions[@]}";
	    do
		if [[ ((${pt} == "65-30" && ${rap} == "0-24"  && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-16"  && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "16-24" && (${bkg} == "pol3") && ${cent} == "0-100") ||
		       (${pt} == "3-30" && ${rap} == "16-24"  && (${bkg} == "pol3") && ${cent} == "0-100") ||
		       (${pt} == "3-65"  && ${rap} == "16-24" && (${bkg} == "pol2") && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-24"  && (${bkg} == "pol1") && ${cent} == "0-20")  ||
		       (${pt} == "65-30" && ${rap} == "0-16"  && (${bkg} == "pol1") && ${cent} == "0-20")  ||
		       (${pt} == "3-30" && ${rap} == "16-24"  && (${bkg} == "pol3") && ${cent} == "0-20")  ||
		       (${pt} == "65-30" && ${rap} == "0-24"  && (${bkg} == "pol0") && ${cent} == "20-40") ||
		       (${pt} == "65-30" && ${rap} == "0-16"  && (${bkg} == "pol1") && ${cent} == "20-40") ||
		       (${pt} == "3-30" && ${rap} == "16-24"  && (${bkg} == "pol2") && ${cent} == "20-40") ||
		       (${pt} == "65-30" && ${rap} == "0-24"  && (${bkg} == "pol0") && ${cent} == "40-100")||
		       (${pt} == "65-30" && ${rap} == "0-16"  && (${bkg} == "pol0") && ${cent} == "40-100")||
		       (${pt} == "3-30" && ${rap} == "16-24"  && (${bkg} == "pol1") && ${cent} == "40-100")
		      ) ]];
		then
		    echo ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt
#		    ratio=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3; printf "\n"; print $5}'`
		    grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3; print $5}'
		    #echo "free alpha"
		    grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeAlpha_allVars.txt | awk '{print $3; print $5}'
		    #echo "free n"
		    grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeN_allVars.txt | awk '{print $3; print $5}'
		    #echo "free nG"
		    grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeGwidth_allVars.txt | awk '{print $3; print $5}'
		    #	    echo ${bkg} ${ratio}
		fi;
	    done;
	done;
    done;
done;

echo "--------------"

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
		if [[ ((${pt} == "65-30" && ${rap} == "0-24"  && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-16"  && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "16-24" && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "3-30"  && ${rap} == "16-24" && (${bkg} == "pol1") && ${cent} == "0-100") ||
		       (${pt} == "3-65"  && ${rap} == "16-24" && (${bkg} == "pol1") && ${cent} == "0-100")
		      ) ]];
		then
		    echo ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt
		   grep 'fracP =' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3; print $5}'
		   if [[ ${bkg} == "pol1" ]];
		   then
		   grep 'fracP =' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeAlpha_allVars.txt | awk '{print $3; print $5}'
		   grep 'fracP =' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeN_allVars.txt | awk '{print $3; print $5}'
		   fi;
		   # echo ${bkg} ${ratio}
		fi;
	    done;
	done;
    done;
done;
