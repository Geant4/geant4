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
#ifndef G4HadronicWhiteBoard_h
#define G4HadronicWhiteBoard_h

#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"

class G4HadronicWhiteBoard
{
  public:
  G4HadronicWhiteBoard(){}
  
  static G4HadronicWhiteBoard & Instance();
  
  void SetProjectile(const G4HadProjectile & aProjectile);
    
  void SetTargetNucleus(const G4Nucleus & aTarget);

  const G4HadProjectile * GetProjectile();
  const G4Nucleus & GetTargetNucleus(); 
  G4ParticleDefinition * GetPDef();
  G4String GetParticleName();
  G4double GetEnergy();
  G4double GetPx();
  G4double GetPy();
  G4double GetPz();
  G4double GetA();
  G4double GetZ();
  
  void Dump();
  
  
  private:
  
  G4HadProjectile * theProjectile;
  G4ParticleDefinition * theDef;
  char * theName;
  G4double theE;
  G4double thePx;
  G4double thePy;
  G4double thePz;
  
  G4Nucleus theTarget;
  G4double theA;
  G4double theZ;
};

#endif
