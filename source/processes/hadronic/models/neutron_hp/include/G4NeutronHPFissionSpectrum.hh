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
// $Id: G4NeutronHPFissionSpectrum.hh,v 1.7 2002-12-12 19:18:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFissionSpectrum_h
#define G4NeutronHPFissionSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPFissionSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPFissionSpectrum()
  {
    expm1 = exp(-1.);
  }
  ~G4NeutronHPFissionSpectrum()
  {
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, eV);
    theThetaDist.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  inline G4double Sample(G4double anEnergy) 
  {
    G4double theta = theThetaDist.GetY(anEnergy);
    // here we need to sample Maxwells distribution, if 
    // need be.
    G4double result, cut;
    G4double range =50*MeV;
    G4double max = Maxwell((theta*eV)/2., theta);
    G4double value;
    do
    {
      result = range*G4UniformRand();
      value = Maxwell(result, theta);
      cut = G4UniformRand();
    }
    while(cut > value/max);
    return result;
  }
  
  private:
 
  // this is the function to sample from. 
  inline G4double Maxwell(G4double anEnergy, G4double theta)
  {
    G4double result = sqrt(anEnergy/eV)*exp(-anEnergy/eV/theta);
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theThetaDist;
  
};

#endif
