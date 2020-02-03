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
/// \file hadronic/Hadr00/include/HistoManager.hh
/// \brief Definition of the HistoManager class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to hold parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Material.hh"

#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4ParticleDefinition;
class HistoManagerMessenger;

class HistoManager
{
public:

  HistoManager();

  ~HistoManager();

  void BeginOfRun();
  void EndOfRun();

  void SetVerbose(G4int val);        

  inline void SetParticleName(const G4String&);
  inline void SetElementName(const G4String&);

  inline void SetNumberOfBinsE(G4int val);
  inline void SetNumberOfBinsP(G4int val);

  inline void SetMinKinEnergy(G4double val);
  inline void SetMaxKinEnergy(G4double val);

  inline void SetMinMomentum(G4double val);
  inline void SetMaxMomentum(G4double val);

  inline void SetHistoName(G4String& val);

  inline void SetTargetMaterial(const G4Material* p);

private:

  HistoManagerMessenger* fMessenger;
  G4AnalysisManager*     fAnalysisManager;

  const G4ParticleDefinition* fNeutron;
  const G4Material* fTargetMaterial;

  G4String fParticleName;
  G4String fElementName;

  G4double fMinKinEnergy;
  G4double fMaxKinEnergy;
  G4double fMinMomentum;
  G4double fMaxMomentum;
 
  G4int fVerbose;
  G4int fBinsE;
  G4int fBinsP;

  G4String fHistoName;
};

inline void HistoManager::SetParticleName(const G4String& name)
{
  fParticleName = name;
}

inline void HistoManager::SetElementName(const G4String& name)
{
  fElementName = name;
}

inline void HistoManager::SetNumberOfBinsE(G4int val)
{
  if(val>0) { fBinsE = val; } 
}

inline void HistoManager::SetNumberOfBinsP(G4int val)
{
  if(val>0) { fBinsP = val; } 
}

inline void HistoManager::SetMinKinEnergy(G4double val)
{
  if(val>0 && val<fMaxKinEnergy) { fMinKinEnergy = val; }
}
 
inline void HistoManager::SetMaxKinEnergy(G4double val)
{
  if(val>fMinKinEnergy) { fMaxKinEnergy = val; }
}

inline void HistoManager::SetMinMomentum(G4double val)
{
  if(val>0 && val<fMaxMomentum) { fMinMomentum = val; }
}
  
inline void HistoManager::SetMaxMomentum(G4double val)
{
  if(val>fMinMomentum) { fMaxMomentum = val; }
}
 
inline void HistoManager::SetHistoName(G4String& val) 
{
  fHistoName = val;
}

inline void HistoManager::SetTargetMaterial(const G4Material* p)
{
  fTargetMaterial = p;
}

#endif
