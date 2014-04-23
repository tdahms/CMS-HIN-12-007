#!/bin/bash

root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_Regit_RD_2013_Trg2_Eff.root\",\"../root_files/Jpsi_Regit_MC_2013_Trg2_Eff.root\",1\)
root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_Regit_RD_2013_Trg3_Eff.root\",\"../root_files/Jpsi_Regit_MC_2013_Trg3_Eff.root\",1\)

root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_Regit_RD_2014_MuId2_Eff.root\",\"../root_files/Jpsi_Regit_MC_2014_MuId2_Eff.root\",1\)
root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_Regit_RD_2014_MuId3_Eff.root\",\"../root_files/Jpsi_Regit_MC_2014_MuId3_Eff.root\",1\)

root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_Regit_RD_2014_Trk_Eff.root\",\"../root_files/Jpsi_Regit_MC_2014_Trk_Eff.root\",1\)



root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_pp_RD_2014_Trg2_Eff.root\",\"../root_files/Jpsi_pp_MC_2014_Trg2_Eff.root\",1\)
root -b -q tnp_trig_eff_correction.C+\(\"../root_files/Jpsi_pp_RD_2014_Trg3_Eff.root\",\"../root_files/Jpsi_pp_MC_2014_Trg3_Eff.root\",1\)

