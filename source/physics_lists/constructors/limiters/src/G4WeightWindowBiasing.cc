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
// ClassName:   G4WeightWindowBiasing
//
// Author:      A. Howard (Nov.22.2013)
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WeightWindowBiasing.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4TransportationManager.hh"
//#include "G4ParallelWorldProcess.hh"
#include "G4WeightWindowProcess.hh"

#include "G4WeightWindowStore.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4WeightWindowBiasing);

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WeightWindowBiasing::G4WeightWindowBiasing(const G4String& name)
  :  G4VPhysicsConstructor(name), fGeomSampler(0), fWWalg(0), fPlaceOfAction(), paraFlag(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WeightWindowBiasing::G4WeightWindowBiasing(G4GeometrySampler* mgs, G4VWeightWindowAlgorithm* wwAlg, G4PlaceOfAction placeOfAction, const G4String& name)
:  G4VPhysicsConstructor(name), fGeomSampler(mgs), fWWalg(wwAlg), fPlaceOfAction(placeOfAction), paraFlag(false), paraName(name)
{
  if(name != "NoParallelWP") {
    paraFlag = true;
    paraName = name;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WeightWindowBiasing::~G4WeightWindowBiasing()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WeightWindowBiasing::ConstructParticle()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WeightWindowBiasing::ConstructProcess()
{
  G4cout << " paraFlag: " << paraFlag << G4endl;

  static G4bool first = true;
  if(first) {
    G4cout << " Preparing WeightWindow Sampling " << G4endl;
    fGeomSampler->SetParallel(paraFlag);
    if(paraFlag) {
      fGeomSampler->PrepareWeightWindow(G4WeightWindowStore::GetInstance(paraName), fWWalg, fPlaceOfAction);
    } else {
      fGeomSampler->PrepareWeightWindow(G4WeightWindowStore::GetInstance(), fWWalg, fPlaceOfAction);
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
