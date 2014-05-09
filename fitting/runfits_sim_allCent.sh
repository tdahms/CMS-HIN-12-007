#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="SimFits_M2242_DblMu0_AllCent_WithSyst_final_forPaper_insert"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

runSyst=0
#ptbins=(6.5-30.0 3.0-30.0 3.0-6.5)
ptbins=(6.5-30.0 3.0-30.0)
#rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)
rapbins=(0.0-1.6 1.6-2.4)
centbins=(0-100 0-20 20-40 40-100)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" ]];
	then
	    bkg="pol1_HI020 pol2_HI2040 pol0_HI40100 pol7";
	elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" ]];
	then
            bkg="pol1_HI020 pol1_HI2040 pol0_HI40100 pol3";
            #bkg="pol1_HI020 pol1_HI2040 pol0_HI40100 pol5";
	    #bkg="pol1_HI020 pol1_HI2040 pol0_HI40100 pol1";
	elif [[ ${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	then
	    bkg="pol3_HI020 pol2_HI2040 pol1_HI40100 pol1";
	    #bkg="pol3_HI020 pol2_HI2040 pol1_HI40100 pol3";
	    #bkg="pol3_HI020 pol2_HI2040 pol2_HI40100 pol2";
	else
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	echo "fitting " ${rap} " " ${pt} " " ${cent} " " ${bkg};
	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 1 &> ${DIRECTORY}/sim_fitCB_polBkg_rap${rap}_pt${pt}_cent${cent}.log
	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 1 &> ${DIRECTORY}/sim_fitCBG_polBkg_rap${rap}_pt${pt}_cent${cent}.log

#	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 1 -a -n -g &> ${DIRECTORY}/sim_fitCB_polBkg_rap${rap}_pt${pt}_cent${cent}_allFree.log
#	./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1CB2 sigCB1CB2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 1 -a -n -g &> ${DIRECTORY}/sim_fitCBG_polBkg_rap${rap}_pt${pt}_cent${cent}_allFree.log

	if [[ ${runSyst} -eq 1 ]];
	then
	    if [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-2.4" ]];
	    then
		bkg="expPol1_HI020 expPol2_HI2040 expPol1_HI40100 expPol3";
	    elif [[ ${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" ]];
	    then
		bkg="expPol1_HI020 expPol1_HI2040 expPol1_HI40100 expPol3";
	    elif [[ ${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		bkg="expPol3_HI020 expPol2_HI2040 expPol1_HI40100 expPol3";
	    else
		echo "skipping " ${rap} " " ${pt};
		continue;
	    fi;
	    echo "fitting " ${rap} " " ${pt} " " ${cent} " " ${bkg};
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 0 &> ${DIRECTORY}/sim_fitCB_expBkg_rap${rap}_pt${pt}_cent${cent}.log
	    ./Fit1DDataPbPbSim -f ../root_files/PbPbData2011_DblMu0_cent0-20_M22-42.root ../root_files/PbPbData2011_DblMu0_cent20-40_M22-42.root ../root_files/PbPbData2011_DblMu0_cent40-100_M22-42.root ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -p ${pt} -y ${rap} -t Mult -l 1 -b 1 -s 0 -u 1 -x 0 -z 0 &> ${DIRECTORY}/sim_fitCBG_expBkg_rap${rap}_pt${pt}_cent${cent}.log
	fi;
    done;
done;

logfiles=(${DIRECTORY}/sim_fitCB*.log)
if [[ -e "${logfiles[0]}" ]];
then
    gzip -f ${DIRECTORY}/sim_fitCB*.log;
fi;
