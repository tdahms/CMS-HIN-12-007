#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"
SUFFIX="M2242_DblMu0_Split"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

runAll=0
#bkgfunctions=(pol0 pol1 pol2 pol3 pol4 pol5 pol6 pol7 expPol1 expPol2 expPol3 expPol4 expPol5 expPol6 expPol7)
bkgfunctions=(pol0 pol1 pol2 pol3 pol4 pol5 pol6 pol7)
#bkgfunctions=(Bern0 Bern1 Bern2 Bern3 Bern4 Bern5 Bern6 Bern7)
#ptbins=(6.5-30.0 3.0-30.0 3.0-6.5)
#rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)
ptbins=(6.5-30.0 3.0-30.0)
rapbins=(0.0-1.6 1.6-2.4)
centbins=(0-100)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${rap} == "0.0-1.6" && ${pt} != "6.5-30.0" ]];
	then
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	if [[ ${rap} == "1.6-2.4" && ${pt} != "3.0-30.0" ]];
	then
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	#	if [[ ${rap} != "1.6-2.4" && ${pt} != "6.5-30.0" ]];
	#	then
	#	    echo "skipping " ${rap} " " ${pt};
	#	    continue;
	#	fi;
	echo "fitting " ${rap} " " ${pt};
	for bkg in "${bkgfunctions[@]}";
	do
	    for bkgP in "${bkgfunctions[@]}";
	    do
		for cent in "${centbins[@]}";
		do
		    if [[ (${runAll} -eq 1 && ${bkg} == ${bkgP}) || (${runAll} -eq 0 &&
		       	  ((${pt} == "6.5-30.0" && ${rap} == "0.0-1.6" && ${bkg} == "pol3" && ${bkgP} == "pol0" && ${cent} == "0-100") ||
			   (${pt} == "3.0-30.0" && ${rap} == "1.6-2.4" && ${bkg} == "pol1" && ${bkgP} == "pol1" && ${cent} == "0-100"))) ]]	
		    then
			echo "fitting " ${cent} " " ${bkg} " " ${bkgP};
			./Fit1DDataPbPbSplit -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v signalCB1 signalCB1P ${bkg} ${bkgP} -d ${DIRECTORY}/ppFracLogCB -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 &> ${DIRECTORY}/pp_fitCB_${bkg}_${bkgP}_rap${rap}_pt${pt}_cent${cent}.log
			./Fit1DDataPbPbSplit -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} ${bkgP} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -b 1 &> ${DIRECTORY}/pp_fitCBG_${bkg}_${bkgP}_rap${rap}_pt${pt}_cent${cent}.log
		    fi;
		done;
	    done;
	done;
    done;
done;

logfiles=(${DIRECTORY}/pp_fitCB*.log)
if [[ -e "${logfiles[0]}" ]];
then
    gzip -f ${DIRECTORY}/pp_fitCB*.log;
fi;

