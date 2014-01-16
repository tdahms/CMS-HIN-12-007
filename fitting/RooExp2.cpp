/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooExp2.cpp,v 1.1 2011/02/01 15:56:28 miheejo Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooExp2 implements a RooResolutionModel that models a exponential
// distribution. Object of class RooExp2 can be used
// for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
// END_HTML
//

#include "TMath.h"

#include "RooFit.h"
#include "Riostream.h"
#include "RooExp2.h"
#include "RooRealConstant.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

//ClassImp(RooExp2);

using namespace std;
using namespace RooFit;

//_____________________________________________________________________________
RooExp2::RooExp2(const char *name, const char *title, RooAbsReal& _xIn,
    RooAbsReal& _cutx, RooAbsReal& _tau1, RooAbsReal& _tau2) :
  RooAbsPdf(name,title), 
  xIn(_xIn.GetName(),_xIn.GetTitle(),this,_xIn),
  cutx("cutx","cutx",this,_cutx),
  tau1("tau1","tau1",this,_tau1),
  tau2("tau2","tau2",this,_tau2)
{  
}

//_____________________________________________________________________________
RooExp2::RooExp2(const RooExp2& other, const char* name) : 
  RooAbsPdf(other,name),
  xIn(other.xIn.GetName(),this,other.xIn),
  cutx("cutx",this,other.cutx),
  tau1("tau1",this,other.tau1),
  tau2("tau2",this,other.tau2)
{
}

//_____________________________________________________________________________
Double_t RooExp2::evaluate() const { 

    Double_t result1 = 0;
    Double_t result2 = 0;
    Double_t otherside = 0;

    if (xIn < cutx) {
      //cout << "x < cutx: " << xIn << " tau1: " << tau1 << endl;
      result1 = exp(tau1*xIn);
    } else {
      //cout << "x >= cutx: " << xIn << " tau2: " << tau2<< endl;
      result1 = exp(tau2*xIn);
      result2 = exp(tau2*cutx);
      otherside = exp(tau1*cutx);
      //cout << "result: " << result1 << endl;
      result1 = result1*(otherside/result2);
    }

    if (xIn < 0.0) result1 = 0.0;

    //cout << "corrected result: " << result1 << endl;
    return result1/fabs(tau1) ; 
}

//_____________________________________________________________________________
Int_t RooExp2::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,xIn)) return 1 ;
  return 0 ;
}


//_____________________________________________________________________________
Double_t RooExp2::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  switch(code) {
  case 1: 
  {
    if (xIn < cutx) {

      Double_t ret(0) ;
      if(tau1 == 0.0) {
        ret = (cutx - xIn.min(rangeName));
      } else {
        ret =  ( exp( tau1*xIn.max(rangeName) ) - exp( tau1*xIn.min(rangeName) ) )/tau1;
      }

      //cout << "Int_exp_dx(c=" << tau1 << ", xmin=" << xIn.min(rangeName) << ", xmax=" << xIn.max(rangeName) << ")=" << ret << endl ;
      return ret ;

    } else {

      Double_t ret(0) ;
      if(tau2 == 0.0) {
        ret = (xIn.max(rangeName) - cutx);
      } else {
        ret =  ( exp( tau2*xIn.max(rangeName) ) - exp( tau2*xIn.min(rangeName) ) )/tau2;
      }

      //cout << "Int_exp_dx(c=" << tau2 << ", xmin=" << xIn.min(rangeName) << ", xmax=" << xIn.max(rangeName) << ")=" << ret << endl ;
      return ret ;

    }
  }
  }
  
  assert(0) ;
  return 0 ;
}
