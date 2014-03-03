#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="SimFits_M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError"

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
ptbins=(6.5-30.0 3.0-30.0)
rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" ]]
	then
	    bkg="pol1_HI020 pol2_HI2040 pol1_HI40100 pol5";
	elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" ]]
	then
	    bkg="pol1_HI020 pol1_HI2040 pol1_HI40100 pol3";
	elif [[ ${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" ]]
	then
	    bkg="pol3_HI020 pol2_HI2040 pol1_HI40100 pol5";
	else
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	echo "fitting " ${rap} " " ${pt};
	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 &> ${DIRECTORY}/sim_fitCB_rap${rap}_pt${pt}.log
	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 &> ${DIRECTORY}/sim_fitCBG_rap${rap}_pt${pt}.log
    done;
done;

logfiles=(${DIRECTORY}/sim_fitCB*.log)
if [[ -e "${logfiles[0]}" ]];
then
    gzip -f ${DIRECTORY}/sim_fitCB*.log;
fi;
