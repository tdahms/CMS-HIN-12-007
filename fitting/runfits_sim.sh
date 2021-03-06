#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="SimFits_M2242_DblMu0_WithSyst_final"

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
ptbins=(6.5-30.0 3.0-6.5)
rapbins=(1.6-2.4)
centbins=(0-100)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	for cent in "${centbins[@]}";
	do
	    if [[ ${cent} != "0-100" && ${pt} != "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		echo "skipping " ${rap} " " ${pt} " " ${cent};
		continue;
	    fi;
	    if [[ ${cent} == "0-100" && ${pt} == "3.0-6.5" && ${rap} == "1.6-2.4" ]];
	    then
		bkg="pol2_HI pol1";
	    elif [[ ${cent} == "0-100" && ${pt} == "6.5-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		bkg="pol3_HI pol1"
	    elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="pol3_HI pol7";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="pol1_HI pol7";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="pol2_HI pol7";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="pol0_HI pol7";
		fi;
	    elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="pol1_HI pol3";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="pol1_HI pol3";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="pol1_HI pol3";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="pol0_HI pol3";
		fi;
	    elif [[ ${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="pol3_HI pol1";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="pol3_HI pol1";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="pol2_HI pol1";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="pol1_HI pol1";
		fi;
	    else
		echo "skipping " ${rap} " " ${pt} " " ${cent};
		continue;
	    fi;
	    echo "fitting " ${rap} " " ${pt} " " ${cent} " " ${bkg};
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 &> ${DIRECTORY}/sim_fitCB_polBkg_rap${rap}_pt${pt}_cent${cent}.log
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 &> ${DIRECTORY}/sim_fitCBG_polBkg_rap${rap}_pt${pt}_cent${cent}.log

	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 -a -n -g &> ${DIRECTORY}/sim_fitCB_polBkg_rap${rap}_pt${pt}_cent${cent}_allFree.log
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 -a -n -g &> ${DIRECTORY}/sim_fitCBG_polBkg_rap${rap}_pt${pt}_cent${cent}_allFree.log

	if [[ ${runSyst} -eq 1 ]];
	then
	    if [[ ${cent} == "0-100" && ${pt} == "3.0-6.5" && ${rap} == "1.6-2.4" ]];
	    then
		bkg="expPol2_HI expPol1";
	    elif [[ ${cent} == "0-100" && ${pt} == "6.5-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		bkg="expPol2_HI expPol1"
	    elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="expPol3_HI expPol3";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="expPol1_HI expPol3";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="expPol2_HI expPol3";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="expPol1_HI expPol3";
		fi;
	    elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="expPol1_HI expPol3";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="expPol1_HI expPol3";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="expPol1_HI expPol3";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="expPol1_HI expPol3";
		fi;
	    elif [[ ${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		if [[ ${cent} == "0-100" ]];
		then
		bkg="expPol2_HI expPol3";
		elif [[ ${cent} == "0-20" ]];
		then
		bkg="expPol3_HI expPol3";
		elif [[ ${cent} == "20-40" ]];
		then
		bkg="expPol2_HI expPol3";
		elif [[ ${cent} == "40-100" ]];
		then
		bkg="expPol1_HI expPol3";
		fi;
	    else
		echo "skipping " ${rap} " " ${pt} " " ${cent};
		continue;
	    fi;
	    echo "fitting " ${rap} " " ${pt} " " ${cent} " " ${bkg};
	    
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 &> ${DIRECTORY}/sim_fitCB_expBkg_rap${rap}_pt${pt}_cent${cent}.log
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t ${cent} -l 1 -b 1 -s 0 -u 1 -x 1 -z 0 &> ${DIRECTORY}/sim_fitCBG_expBkg_rap${rap}_pt${pt}_cent${cent}.log
	fi;
	done;
    done;
done;

logfiles=(${DIRECTORY}/sim_fitCB*.log)
if [[ -e "${logfiles[0]}" ]];
then
    gzip -f ${DIRECTORY}/sim_fitCB*.log;
fi;
