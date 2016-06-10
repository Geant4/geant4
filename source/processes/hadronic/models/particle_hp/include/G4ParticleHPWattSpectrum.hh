//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPWattSpectrum_h
#define G4ParticleHPWattSpectrum_h 1

#include <fstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4ParticleHPVector.hh"
#include "G4VParticleHPEDis.hh"

// we will need a List of these .... one per term.

class G4ParticleHPWattSpectrum : public G4VParticleHPEDis
{
  public:
  G4ParticleHPWattSpectrum()
  {
    expm1 = G4Exp(-1.);
  }
  ~G4ParticleHPWattSpectrum()
  {
  }
  
  inline void Init(std::istream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, CLHEP::eV);
    theApar.Init(aDataFile, CLHEP::eV);
    theBpar.Init(aDataFile, CLHEP::eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy);
  
  private:
  
  inline G4double Watt(G4double anEnergy, G4double a, G4double b)
  {
    G4double energy = anEnergy/CLHEP::eV;
    G4double result = G4Exp(-energy/a)*std::sinh(std::sqrt(b*energy));
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4ParticleHPVector theFractionalProb;
  
  G4ParticleHPVector theApar;
  G4ParticleHPVector theBpar;
  
};

#endif
