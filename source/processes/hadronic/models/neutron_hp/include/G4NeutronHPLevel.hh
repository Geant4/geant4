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
// $Id: G4NeutronHPLevel.hh,v 1.5 2001-07-11 10:07:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLevel_h
#define G4NeutronHPLevel_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4DynamicParticleVector.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
class G4NeutronHPGamma;

class G4NeutronHPLevel
{
  public:
  
  G4NeutronHPLevel() 
  {
    nGammas = 0;
    theGammas = NULL;
  }

  ~G4NeutronHPLevel();
  
  void SetNumberOfGammas(G4int aGammas);
  
  void SetGamma(G4int i, G4NeutronHPGamma * aGamma);
  
  G4DynamicParticleVector * GetDecayGammas();
    
  inline void SetLevelEnergy(G4double anEnergy)
  {
    levelEnergy = anEnergy;
  }
  
  inline G4double GetLevelEnergy()
  {
    return levelEnergy;
  }

  G4double GetGammaEnergy(G4int i);
  
  private:
  
  G4double levelEnergy;  

  G4int nGammas;
  G4NeutronHPGamma ** theGammas;
};

#endif
