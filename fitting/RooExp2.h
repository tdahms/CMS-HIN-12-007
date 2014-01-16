/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooExp2.h,v 1.1 2011/02/01 15:56:33 miheejo Exp $
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
#ifndef ROO_Exp2
#define ROO_Exp2

#include "Riostream.h"

#include <math.h>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooExp2 : public RooAbsPdf {
public:

  // Constructors, assignment etc
  RooExp2() { }
  RooExp2(const char *name, const char *title, RooAbsReal& x, 
		RooAbsReal& cutx, RooAbsReal& tau1, RooAbsReal& tau2 );

  RooExp2(const RooExp2& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooExp2(*this,newname) ; }
  inline virtual ~RooExp2() {}
  
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;

protected:

  Double_t evaluate() const ;
//  virtual Double_t evaluate() const ;

  RooRealProxy xIn ;
  RooRealProxy cutx ;
  RooRealProxy tau1 ;
  RooRealProxy tau2 ;

private:
//  ClassDef(RooExp2,1)
};

#endif
