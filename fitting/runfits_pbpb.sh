#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="DblMu0"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-1.6 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-1.6 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent0-20.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 0-20 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-1.6 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-1.6 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent20-40.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 20-40 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 0.0-1.6 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 0.0-1.6 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v signalCB1 signalCB1P pol1  -d ${DIRECTORY}/fracLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/PbPbData2011_DblMu0_cent40-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/fracLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 40-100 -x 4 -r 1 -l 1 -s 0 -n 1
