#!/bin/bash

#root -b -q make_psi2s_mupt.C+\(6.5,30.0,0.0,2.4,1,1\)
#root -b -q make_psi2s_mupt.C+\(6.5,30.0,0.0,1.6,1,1\)
#root -b -q make_psi2s_mupt.C+\(6.5,30.0,1.6,2.4,1,1\)
#root -b -q make_psi2s_mupt.C+\(3.0,30.0,1.6,2.4,1,1\)
#root -b -q make_psi2s_mupt.C+\(3.0,6.5,1.6,2.4,1,1\)

root -b -q make_Jpsi_mupt.C+\(3.0,30.0,1.6,2.4,1,1\)
root -b -q make_Jpsi_mupt.C+\(6.5,30.0,0.0,2.4,1,1\)
root -b -q make_Jpsi_mupt.C+\(6.5,30.0,0.0,1.6,1,1\)
root -b -q make_Jpsi_mupt.C+\(6.5,30.0,1.6,2.4,1,1\)
root -b -q make_Jpsi_mupt.C+\(3.0,6.5,1.6,2.4,1,1\)
