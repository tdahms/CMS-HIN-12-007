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
		# exclude fits gone bad
		if [[ (${pt} == "65-30" && ${rap} == "0-24"  && ${bkg} == "gausBkg" && ${cent} == "0-100") ||
		      (${pt} == "3-65"  && ${rap} == "16-24" && ${bkg} == "pol6"    && ${cent} == "0-100") ||
		      (${pt} == "65-30" && ${rap} == "0-16"  && ${bkg} == "gausBkg" && ${cent} == "0-20") ||
		      (${bkg} == "pol6") || (${bkg} == "pol5") ]]
		then
#		    echo "skipping rap" ${rap} " pT"${pt} " cent" ${cent} " " ${bkg};
		    continue;
		fi;
#		echo "reading rap" ${rap} " pT"${pt} " cent" ${cent} " " ${bkg};
		ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
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
		    ratio=${ratio_tmp}
		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeAlpha_allVars.txt | awk '{print $3}'`
		    echo "freeAlpha" ${ratio_tmp} >> fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt

		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeN_allVars.txt | awk '{print $3}'`
		    echo "freeN" ${ratio_tmp} >> fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt

		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeGwidth_allVars.txt | awk '{print $3}'`
		    echo "freeGwidth" ${ratio_tmp} >> fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
#		    echo "CB" ${ratio_tmp} >> fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
		fi;
	    done;
	    sort -k2 fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt > fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt
	    rm -f fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
	    minBkg=`head -n1 fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $1}'`
	    maxBkg=`tail -n1 fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $1}'`
	    min=`head -n1 fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $2}'`
	    max=`tail -n1 fracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $2}'`
	    spread=`echo "(${max}-(${min}))*100"|bc`
	    hi=`echo "(${max}-(${ratio}))*100"|bc`
	    lo=`echo "(${ratio}-(${min}))*100"|bc`
	    relRatio=`echo "scale=1;${spread}/${ratio}"|bc`
	    relRatioHi=`echo "scale=1;${hi}/${ratio}"|bc`
	    relRatioLo=`echo "scale=1;${lo}/${ratio}"|bc`
	    echo "PbPb rap" ${rap} " pT"${pt} " cent" ${cent}": " ${relRatio}"% --- +"${relRatioHi}"/-"${relRatioLo}"% ("${maxBkg} ${minBkg}")"
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
		# exclude fits gone bad
		if [[ (${pt} == "3-65" && ${rap} == "16-24" && ${bkg} == "gausBkg" && ${cent} == "0-100") ]]
		then
		    continue;
		fi;
		ratio_tmp=`grep 'fracP =' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
		echo ${bkg} ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
		if [[ ((${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "0-16" && ${bkg} == "pol3" && ${cent} == "0-100") ||
		       (${pt} == "65-30" && ${rap} == "16-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3-30"  && ${rap} == "16-24" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3-65"  && ${rap} == "16-24" && ${bkg} == "pol1" && ${cent} == "0-100")
		      ) ]];
		then
		    ratio=${ratio_tmp}
		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeAlpha_allVars.txt | awk '{print $3}'`
		    echo "freeAlpha" ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt

		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeN_allVars.txt | awk '{print $3}'`
		    echo "freeN" ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt

		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_freeGwidth_allVars.txt | awk '{print $3}'`
		    echo "freeGwidth" ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt

#		    ratio_tmp=`grep 'fracP = ' ${DIRECTORY}/ppFracLogCB_${bkg}_rap${rap}_pT${pt}_cent${cent}_allVars.txt | awk '{print $3}'`
#		    echo "CB" ${ratio_tmp} >> ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
		fi;
	    done;
	    sort -k2 ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt > ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt
	    rm -f ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP.txt
	    minBkg=`head -n1 ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $1}'`
	    maxBkg=`tail -n1 ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $1}'`
	    min=`head -n1 ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $2}'`
	    max=`tail -n1 ppFracLogCBG_rap${rap}_pT${pt}_cent${cent}_fracP_sorted.txt|awk '{print $2}'`
	    spread=`echo "(${max}-(${min}))*100"|bc`
	    hi=`echo "(${max}-(${ratio}))*100"|bc`
	    lo=`echo "(${ratio}-(${min}))*100"|bc`
	    relRatio=`echo "scale=1;${spread}/${ratio}"|bc`
	    relRatioHi=`echo "scale=1;${hi}/${ratio}"|bc`
	    relRatioLo=`echo "scale=1;${lo}/${ratio}"|bc`
	    echo "pp rap" ${rap} " pT"${pt} " cent" ${cent}": " ${relRatio}"% --- +"${relRatioHi}"/-"${relRatioLo}"% ("${maxBkg} ${minBkg}")"
	done;
    done;
done;
