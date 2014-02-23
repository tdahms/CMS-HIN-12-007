#! /bin/bash


DATE="20140221"
SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError"

DIRECTORY="${DATE}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "${DIRECTORY} does not exist."
    exit 1
else
    echo "Reading results in ${DIRECTORY}"
fi

bkgfunctions=(pol1 pol2 pol3 pol4 pol5 pol6 expFunct gausBkg)
ptbins=(65-30 3-30 3-65)
rapbins=(0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)

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
		if [[ (${pt} == "65-30" && ${rap} == "0-24" && ${bkg} == "pol1" && ${cent} == "0-100") ||
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
		       ]];
		then
		    echo "Running C&F scan for ${pt} ${rap} ${cent} ${bkg}";
		    root -b -q StandardHypoTestInvDemo.C+\(\"${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_Workspace.root\",\"workspace\",\"model\",\"\",\"data\",0,2,false,30,0,0.15,10000,false\)
		fi;
	    done;
	done;
    done;
done;
