#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

runSyst=1
bkgfunctions=(pol1 pol2 pol3 pol4 pol5 pol6 expFunct gausBkg)
ptbins=(6.5-30.0 3.0-30.0 3.0-6.5)
rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)
centbins=(0-100)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${rap} != "1.6-2.4" && ${pt} != "6.5-30.0" ]];
	then
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	echo "fitting " ${rap} " " ${pt};
	for bkg in "${bkgfunctions[@]}";
	do
	    for cent in "${centbins[@]}";
	    do
		echo "fitting " ${cent} " " ${bkg};
		./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 &> ${DIRECTORY}/pp_fitCB_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
		./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
		if [[ ${runSyst} -eq 1 &&
		      ((${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" && ${bkg} == "pol3" && ${cent} == "0-100") ||
		       (${pt} == "6.5-30.0" && ${rap} == "1.6-2.4" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" && ${bkg} == "pol5" && ${cent} == "0-100") ||
		       (${pt} == "3.0-6.5"  && ${rap} == "1.6-2.4" && ${bkg} == "pol1" && ${cent} == "0-100")
		      ) ]];
		then
		    echo "Running systematic fit for ${pt} ${rap} ${bkg}";
		    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
		    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
		    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log
		fi;
	    done;
	done;
    done;
done;

logfiles=(${DIRECTORY}/pp_fitCB*.log)
if [[ -e "${logfiles[0]}" ]];
then
    gzip -f ${DIRECTORY}/pp_fitCB*.log;
fi;

