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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4SCProcessorStand
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 10.05.2002
//
// Modifications: 
//
// 09.12.2002 remove warning (VI)
// 23.12.2002 change interface in order to move to cut per region (VI)
// 26-12-02 Secondary production moved to derived classes (VI)
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4SCProcessorStand.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4EmModelManager.hh"
#include "G4MaterialTable.hh"
#include "G4VEmModel.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SCProcessorStand::G4SCProcessorStand() 
  : G4VSubCutoffProcessor(),
    theLambdaSubTable(0),
    thePositron(G4Positron::Positron()),
    modelManager(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4SCProcessorStand::~G4SCProcessorStand() 
{
  if(theLambdaSubTable) theLambdaSubTable->clearAndDestroy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4SCProcessorStand::Initialise(const G4ParticleDefinition* p,
                                    const G4ParticleDefinition* sp, 
                                          G4EmModelManager* m) 
{
  particle = p;
  secondaryParticle = sp;
  modelManager = m;
  navigator = (G4TransportationManager::GetTransportationManager())
                                       ->GetNavigatorForTracking();
  theCuts = modelManager->Cuts();
  theSubCuts = modelManager->Cuts();
  initialMass= particle->GetPDGMass();

  // Access to materials
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nMaterials = G4Material::GetNumberOfMaterials();

  rangeCuts.clear();

  for(size_t i=0; i<nMaterials; i++) {

    const G4Material* aMaterial = (*theMaterialTable)[i];
    rangeCuts.push_back(secondaryParticle->GetRangeThreshold(aMaterial));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4std::vector<G4Track*>*  G4SCProcessorStand::SampleSecondaries(
                    const G4Step& step,
                    const G4DynamicParticle* dp,
			  G4double tmax,
                          G4double meanLoss)
{
  G4bool b;
  const G4Material* mat = (step.GetTrack())->GetMaterial();

  if(mat != material) {
    material = mat; 
    materialIndex = material->GetIndex();
    cut    = (*theCuts)[materialIndex];
    subcut = (*theSubCuts)[materialIndex];
    rcut   = rangeCuts[materialIndex];
  }

  if(subcut >= cut) return 0;
  tmax = G4std::min(tmax,subcut);

  G4double ekin  = dp->GetKineticEnergy();
  G4double effChargeFactor = 1.0;
  G4double massRatio = 1.0;
  G4double mass = initialMass;

  if(dp->GetDefinition() != particle) {
    mass = dp->GetMass();
    massRatio = initialMass/mass;
    G4double q = particle->GetPDGCharge()/dp->GetCharge();
    effChargeFactor = q*q;
  } 

  G4double cross = (*theLambdaSubTable)[materialIndex]->GetValue(ekin*massRatio, b);

  if(0.0 >= cross) return 0;

  G4StepPoint* pre = step.GetPreStepPoint();
  G4double presafety  = pre->GetSafety();
  G4ThreeVector postpoint = step.GetPostStepPoint()->GetPosition();
  G4double postsafety = navigator->ComputeSafety(postpoint);
  G4double safety = G4std::min(presafety,postsafety);
  if(safety >= rcut) return 0;
 

  G4ThreeVector prepoint = pre->GetPosition();
  G4ThreeVector dr = postpoint - prepoint;
  G4double pretime = step.GetPreStepPoint()->GetGlobalTime();
  G4double fragment = 0.0;
  G4double dt = 0.0;
  G4double length = step.GetStepLength();
  G4double inv_v = (ekin + mass)/(c_light*dp->GetTotalMomentum());

  G4VEmModel* currentModel = modelManager->SelectModel(ekin);
  G4std::vector<G4Track*>* vtr = new G4std::vector<G4Track*>;

  do {

    G4double del = G4UniformRand()*effChargeFactor / cross;
    fragment += del/length;
    if (fragment > 1.0) break;
             
    dt += del * inv_v; 
    G4std::vector<G4DynamicParticle*>* newp = 
           currentModel->SampleSecondaries(material, dp, subcut, cut);
    if (newp) {

      G4DynamicParticle* p;
      G4int nNew = newp->size();
      for (G4int i=0; i<nNew; i++) {

	p = (*newp)[i];
        G4double e = p->GetKineticEnergy();
        if (p->GetDefinition() == thePositron) e += electron_mass_c2;

        if (e <= meanLoss) {
          meanLoss  -= e;
          G4Track* t = new G4Track(p, pretime + dt, prepoint + fragment*dr);
          vtr->push_back(t);

        } else {

          delete p;
        }
      }
    }
  } while (fragment < 1.0);

  return vtr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 




