//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst33Sigma.hh,v 1.1 2002-10-29 15:43:05 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33Sigma
//
// Class description:
//
// A class for simple Gaussian statistics taking into account 
// weights.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef Tst33Sigma_hh
#define Tst33Sigma_hh Tst33Sigma_hh

#include "g4std/iostream"
#include "globals.hh"

class Tst33Sigma
{ 
  
 public:  // with description

  Tst33Sigma();
    // call Init()

  ~Tst33Sigma();
    // simple deonstruction

  void Init();
    // initialisation

  G4int Xin(G4double x, G4double w = 1.);
   // enter a value [with a weight] 

  G4int GetEntries() const;
    // get the number of entries

  G4double GetMean() const;
    // get the weighted mean value: Sum(w*x)/Sum(w)

  G4double GetSigma() const;
    // get the weighted sigma: sqrt(Sum(w*x*x)/Sum(w)-mean^2)

  G4double GetXsum() const;
    // get sum over x: Sum(x) 

  G4double GetXXsum() const;
    // get sum over squared x: Sum(x*x)

  G4double GetSumOfWeights() const;
    // get sum over weights: Sum(w) 

  G4double GetWeightedXsum() const;
    // get weighted sum over x*w: Sum(w*x)

  G4double GetWeightedXXsum() const;
    // get weighted sum over squarted x: Sum(W*x*x)

 private:

  G4int GetCalc() const;
  G4int Calculate() const;
  void Error(const G4String &m);

 private:

  G4int fEntries;
  mutable G4double fMean;
  mutable G4double fSigma;
  G4double fXsum;
  G4double fXXsum;
  G4double fWsum;
  G4double fWXsum;
  G4double fWXXsum;
  mutable G4int fcalc;

};

G4std::ostream& operator<<(G4std::ostream &out, const Tst33Sigma &s);

#endif
