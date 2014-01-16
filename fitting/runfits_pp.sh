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

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol1 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol1 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol2 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol2 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol2 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol2 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol2 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol2 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol2 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol2 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol2 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol2 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol3 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol3 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol3 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol3 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol3 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol3 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol3 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol3 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol3 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol3 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol4 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol4 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol4 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol4 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol4 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol4 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol4 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol4 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol4 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol4 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol5 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol5 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol5 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol5 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol5 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol5 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol5 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol5 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol5 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol5 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1

./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol6 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol6 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol6 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol6 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 0.0-1.6 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol6 -d ${DIRECTORY}/ppFracLogCB -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol6 -d ${DIRECTORY}/ppFracLogCBG -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol6 -d ${DIRECTORY}/ppFracLogCB -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol6 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v signalCB1 signalCB1P pol6 -d ${DIRECTORY}/ppFracLogCB -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent0-100.root -v sigCB1G2 sigCB1G2P pol6 -d ${DIRECTORY}/ppFracLogCBG -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4 -r 1 -l 1 -s 0 -n 1
