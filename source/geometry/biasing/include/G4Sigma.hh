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
// $Id: G4Sigma.hh,v 1.5 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Sigma
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4Sigma_hh
#define G4Sigma_hh

#include "g4std/iostream"
#include "globals.hh"

class G4Sigma
{ 
  
 public:  // with description

  G4Sigma();
  ~G4Sigma();

  void Init();
  G4int Xin(G4double x, G4double w = 1.);
  G4int GetEntries() const;
  G4double GetMean() const;
  G4double GetSigma() const;
  G4double GetXsum() const;
  G4double GetXXsum() const;
  G4double GetSumOfWeights() const;
  G4double GetWeightedXsum() const;
  G4double GetWeightedXXsum() const;

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

G4std::ostream& operator<<(G4std::ostream &out, const G4Sigma &s);

#endif
