#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

bkgfunctions=(pol1 pol2 pol3 pol4 pol5 pol6)
ptbins=(6.5-30.0 3.0-30.0 3.0-6.5)
rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)
#rapbins=(1.6-2.4)
centbins=(0-100 0-20 20-40 40-100)

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
	for cent in "${centbins[@]}";
	do
	    if [[ ${cent} != "0-100" && ${pt} != "3.0-30.0" && ${rap} == "1.6-2.4" ]];
	    then
		echo "skipping " ${rap} " " ${pt} " " ${cent};
		continue;
	    fi;
	    for bkg in "${bkgfunctions[@]}";
	    do
		echo "fitting " ${cent} " " ${bkg};
		./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/fracLogCB  -u -p ${pt} -y ${rap} -t ${cent} -x 4 -r 1 -l 1 -s 0 -n 1 &> ${DIRECTORY}/pbpb_fitCB_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
		./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/fracLogCBG -u -p ${pt} -y ${rap} -t ${cent} -x 4 -r 1 -l 1 -s 0 -n 1 &> ${DIRECTORY}/pbpb_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
	    done;
	done;
    done;
done;

gzip -f ${DIRECTORY}/pbpb_fitCB*.log
