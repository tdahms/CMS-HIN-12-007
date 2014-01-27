#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="M2242_DblMu0"

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

for bkg in "${bkgfunctions[@]}";
do
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
    ./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
done;
