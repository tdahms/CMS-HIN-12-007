#! /bin/bash

#low pt
#0-100%
N=1000
parVals="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap1624_pt330_cent0-100_output_0.txt"
randfname="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root"
systfname2="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root"
systfname3="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_preFitResult.root"
systErr=0.33
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#0-20%
N=1000
parVals="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-20_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap1624_pt330_cent0-20_output_0.txt"
randfname="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-20_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root" # fix alpha and n to 0-100
systfname2="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-20_dPhi_fitResult.root" # fix cutx to 0-20
systfname3="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent0-20_dPhi_preFitResult.root" # fix gaus to prefit 0-20
systErr=0.33
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#20-40%
N=100
parVals="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent20-40_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap1624_pt330_cent20-40_output.txt"
randfname="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent20-40_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root" # fix alpha and n to 0-100
systfname2="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent20-40_dPhi_fitResult.root" # fix cutx to 20-40
systfname3="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent20-40_dPhi_preFitResult.root" # fix gaus to prefit 20-40
systErr=0.32
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#40-100%
N=100
parVals="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent40-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap1624_pt330_cent40-100_output.txt"
randfname="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent40-100_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent0-100_dPhi_fitResult.root" # fix alpha and n to 0-100
systfname2="20120524_simFits/simLogCB_rap1.6-2.4_pT3.0-30.0_cent40-100_dPhi_fitResult.root" # fix cutx to 40-100
systfname3="20120524_simFits/simLogCBG_rap1.6-2.4_pT3.0-30.0_cent40-100_dPhi_preFitResult.root" # fix gaus to prefit 40-100
systErr=0.34
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)

#high pt
#0-100%
N=1000
parVals="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap024_pt6530_cent0-100_output_0.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname2="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname3="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-100_dPhi_preFitResult.root"
systErr=0.14
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#0-20%
N=1000
parVals="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-20_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap024_pt6530_cent0-20_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-20_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname2="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent0-20_dPhi_fitResult.root"
systfname3="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent0-20_dPhi_preFitResult.root"
systErr=0.10
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#20-40%
N=100
parVals="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent20-40_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap024_pt6530_cent20-40_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent20-40_dPhi_fitResult.root"
systfname="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent20-40_dPhi_fitResult.root"
systErr=0.10
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#40-100%
N=100
parVals="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent40-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap024_pt6530_cent40-100_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-2.4_pT6.5-30.0_cent40-100_dPhi_fitResult.root"
systfname="20120524_simFits/simLogCB_rap0.0-2.4_pT6.5-30.0_cent40-100_dPhi_fitResult.root"
systErr=0.15
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#high pt mid-rapidity
#0-100%
N=1000
parVals="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap016_pt6530_cent0-100_output_0.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname2="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname3="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-100_dPhi_preFitResult.root"
systErr=0.14
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#0-20%
N=1000
parVals="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-20_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap016_pt6530_cent0-20_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-20_dPhi_fitResult.root"
systfname1="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent0-100_dPhi_fitResult.root"
systfname2="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent0-20_dPhi_fitResult.root"
systfname3="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent0-20_dPhi_preFitResult.root"
systErr=0.10
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#20-40%
N=100
parVals="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent20-40_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap016_pt6530_cent20-40_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent20-40_dPhi_fitResult.root"
systfname="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent20-40_dPhi_fitResult.root"
systErr=0.10
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)


#40-100%
N=100
parVals="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent40-100_FitResult.txt"
outfname="20120524_simFits/toyMC_simFit_rap016_pt6530_cent40-100_output.txt"
randfname="20120524_simFits/simLogCBG_rap0.0-1.6_pT6.5-30.0_cent40-100_dPhi_fitResult.root"
systfname="20120524_simFits/simLogCB_rap0.0-1.6_pT6.5-30.0_cent40-100_dPhi_fitResult.root"
systErr=0.15
seed=0

root -b -q  toyMC.cpp+\($N,\"$parVals\",\"$outfname\",\"$randfname\",\"$systfname1\",\"$systfname2\",\"$systfname3\",1.0,$seed,1.0,$systErr\)
