#!/bin/bash

root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,0.0,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,0.0,1.6,1,1,1\)
root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,1.6,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(1,3.0,30.0,1.6,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(1,3.0,6.5,1.6,2.4,1,1,1\)

root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,0.0,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,0.0,1.6,1,1,1\)
root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,1.6,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(0,3.0,30.0,1.6,2.4,1,1,1\)
root -b -q .x fitMC_CBG.C+\(0,3.0,6.5,1.6,2.4,1,1,1\)

root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,0.0,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,0.0,1.6,1,0,1\)
root -b -q .x fitMC_CBG.C+\(1,6.5,30.0,1.6,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(1,3.0,30.0,1.6,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(1,3.0,6.5,1.6,2.4,1,0,1\)

root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,0.0,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,0.0,1.6,1,0,1\)
root -b -q .x fitMC_CBG.C+\(0,6.5,30.0,1.6,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(0,3.0,30.0,1.6,2.4,1,0,1\)
root -b -q .x fitMC_CBG.C+\(0,3.0,6.5,1.6,2.4,1,0,1\)
