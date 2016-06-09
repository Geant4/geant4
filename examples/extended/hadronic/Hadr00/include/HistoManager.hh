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
// $Id: HistoManager.hh,v 1.2 2010-10-11 11:02:36 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Histo;
class G4ParticleDefinition;

class HistoManager
{
public:

  static HistoManager* GetPointer();

private:

  HistoManager();

public: 

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

  inline G4int GetVerbose() const;

private:

  static HistoManager* fManager;

  Histo* histo;

  const G4ParticleDefinition* neutron;

  G4String particleName;
  G4String elementName;

  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double minMomentum;
  G4double maxMomentum;
 
  G4int verbose;
  G4int nBinsE;
  G4int nBinsP;

  G4bool isInitialised;
};

inline void HistoManager::SetParticleName(const G4String& name)
{
  particleName = name;
}

inline void HistoManager::SetElementName(const G4String& name)
{
  elementName = name;
}

inline void HistoManager::SetNumberOfBinsE(G4int val)
{
  if(val>0) { nBinsE = val; } 
}

inline void HistoManager::SetNumberOfBinsP(G4int val)
{
  if(val>0) { nBinsP = val; } 
}

inline void HistoManager::SetMinKinEnergy(G4double val)
{
  if(val>0 && val<maxKinEnergy) { minKinEnergy = val; }
}
 
inline void HistoManager::SetMaxKinEnergy(G4double val)
{
  if(val>minKinEnergy) { maxKinEnergy = val; }
}

inline void HistoManager::SetMinMomentum(G4double val)
{
  if(val>0 && val<maxMomentum) { minMomentum = val; }
}
  
inline void HistoManager::SetMaxMomentum(G4double val)
{
  if(val>minMomentum) { maxMomentum = val; }
}
 
inline G4int HistoManager::GetVerbose() const
{
  return verbose;
}

#endif
