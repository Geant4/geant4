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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4ImportanceBiasing
//
// Author:      A. Howard (Nov.09.2013)
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ImportanceBiasing.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4TransportationManager.hh"
//#include "G4ParallelWorldProcess.hh"
#include "G4ImportanceProcess.hh"

#include "G4IStore.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4ImportanceBiasing);

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ImportanceBiasing::G4ImportanceBiasing(const G4String& name)
:  G4VPhysicsConstructor(name), fGeomSampler(0), paraFlag(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ImportanceBiasing::G4ImportanceBiasing(G4GeometrySampler* mgs, const G4String& name)
:  G4VPhysicsConstructor(name), fGeomSampler(mgs), paraFlag(false), paraName(name)
{
  if(name != "NoParallelWP") {
    paraFlag = true;
    paraName = name;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ImportanceBiasing::~G4ImportanceBiasing()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ImportanceBiasing::ConstructParticle()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ImportanceBiasing::ConstructProcess()
{
  G4cout << " paraFlag: " << paraFlag << G4endl;

  static G4bool first = true;
  if(first) {
    G4cout << " Preparing Importance Sampling " << G4endl;
    fGeomSampler->SetParallel(paraFlag);
    if(paraFlag) {
      fGeomSampler->PrepareImportanceSampling(G4IStore::GetInstance(paraName), 0);
    } else {
      fGeomSampler->PrepareImportanceSampling(G4IStore::GetInstance(), 0);
    }
  }

  if(first) {
    fGeomSampler->Configure();
    first = false;
  }

#ifdef G4MULTITHREADED
  fGeomSampler->AddProcess();
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
