#! /bin/bash

root -b -q fitSideBands.C+\(0,3.0,30.0,1.6,2.4,1,1,0,1,1.5\)
root -b -q fitSideBands.C+\(1,3.0,30.0,1.6,2.4,1,1,0,3,1.5\)

root -b -q fitSideBands.C+\(0,6.5,30.0,0.0,1.6,1,1,0,3,1.5\)
root -b -q fitSideBands.C+\(1,6.5,30.0,0.0,1.6,1,1,0,1,1.5\)
