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
#ifndef G4HadronicWhiteBoard_h
#define G4HadronicWhiteBoard_h

#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"


class G4HadronicWhiteBoard
{
  public:
  G4HadronicWhiteBoard();
  
  static G4HadronicWhiteBoard& Instance();
  
  void SetProjectile(const G4HadProjectile& aProjectile);
    
  void SetTargetNucleus(const G4Nucleus& aTarget);

  void SetProcessName(const G4String& aProcessName);

  void SetModelName(const G4String& aModelName);

  const G4HadProjectile* GetProjectile();
  const G4Nucleus& GetTargetNucleus(); 
  const G4ParticleDefinition* GetPDef();
  G4String GetParticleName();
  G4double GetEnergy();
  G4double GetPx();
  G4double GetPy();
  G4double GetPz();
  G4int GetA();
  G4int GetZ();

  void Dump();
  
  
  private:

  static G4ThreadLocal G4HadronicWhiteBoard *theInstance;

  const G4HadProjectile* theProjectile;
  const G4ParticleDefinition* theDef;
  const char* theName;
  G4double theE;
  G4double thePx;
  G4double thePy;
  G4double thePz;
  
  G4Nucleus theTarget;
  G4int theA;
  G4int theZ;

  G4String theProcessName;
  G4String theModelName;
};

#endif
