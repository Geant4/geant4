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
// $Id: G4NeutronHPWattSpectrum.hh,v 1.7 2002-12-12 19:18:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPWattSpectrum_h
#define G4NeutronHPWattSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPWattSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPWattSpectrum()
  {
    expm1 = exp(-1.);
  }
  ~G4NeutronHPWattSpectrum()
  {
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, eV);
    theApar.Init(aDataFile, eV);
    theBpar.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy);
  
  private:
  
  inline G4double Watt(G4double anEnergy, G4double a, G4double b)
  {
    G4double energy = anEnergy/eV;
    G4double result = exp(-energy/a)*sinh(sqrt(b*energy));
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theApar;
  G4NeutronHPVector theBpar;
  
};

#endif
