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
// $Id: G4WentzelVIRelModel.cc 104487 2017-06-01 13:41:55Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelVIRelModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 08.06.2012 from G4WentzelVIRelModel
//
// Modifications:
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// G.Wentzel, Z. Phys. 40 (1927) 590.
// H.W.Lewis, Phys Rev 78 (1950) 526.
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
// L.Urban, CERN-OPEN-2006-077.

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIRelModel.hh"
#include "G4WentzelVIRelXSection.hh"
#include "G4WentzelOKandVIxSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4double> G4WentzelVIRelModel::effMass;

#ifdef G4MULTITHREADED
G4Mutex G4WentzelVIRelModel::WentzelVIRelModelMutex;
#endif

G4WentzelVIRelModel::G4WentzelVIRelModel() :
  G4WentzelVIModel(true, "WentzelVIRel")
{
  fNistManager = G4NistManager::Instance();
  G4WentzelVIRelXSection* ptr = new G4WentzelVIRelXSection();
  SetWVICrossSection(static_cast<G4WentzelOKandVIxSection*>(ptr));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIRelModel::~G4WentzelVIRelModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIRelModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector& cuts)
{
  // Access to materials
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  if(theCoupleTable->GetTableSize() != effMass.size()) {
    ComputeEffectiveMass();
  }

  G4WentzelVIModel::Initialise(p, cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4WentzelVIRelModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    SetCurrentCouple(cup); 
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
    GetWVICrossSection()->SetTargetMass(effMass[currentMaterialIndex]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
                             G4double kinEnergy,
                             G4double Z, G4double,
                             G4double cutEnergy, G4double)
{
  G4double cross = 0.0;
  if(p != particle) { SetupParticle(p); }
  if(kinEnergy < lowEnergyLimit) { return cross; }
  if(!CurrentCouple()) {
    G4Exception("G4WentzelVIRelModel::ComputeCrossSectionPerAtom", "em0011",
                FatalException, " G4MaterialCutsCouple is not defined");
    return cross;
  }
  DefineMaterial(CurrentCouple());
  G4int iz = G4lrint(Z);
  G4double tmass = (1 == iz) ? CLHEP::proton_mass_c2
    : fNistManager->GetAtomicMassAmu(iz)*amu_c2;
  wokvi->SetTargetMass(tmass);
  cosTetMaxNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);
  if(cosTetMaxNuc < 1.0) {
    G4double cost = wokvi->SetupTarget(iz, cutEnergy);
    cross = wokvi->ComputeTransportCrossSectionPerAtom(cost);
    /*    
    //if(p->GetParticleName() == "e-")      
    G4cout << "G4WentzelVIRelModel::CS: Z= " << G4int(Z) 
           << " e(MeV)= " << kinEnergy 
           << " 1-cosN= " << 1 - cosTetMaxNuc << " cross(bn)= " << cross/barn
           << " " << particle->GetParticleName() << G4endl;
    */
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIRelModel::ComputeEffectiveMass()
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4WentzelVIRelModel::WentzelVIRelModelMutex);
#endif
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t ncouples = theCoupleTable->GetTableSize();
  if(ncouples != effMass.size()) {
    effMass.resize(ncouples, 0.0);
    for(size_t i=0; i<ncouples; ++i) {
      const G4Material* mat = 
	theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
      const G4ElementVector* elmVector = mat->GetElementVector();
      G4int nelm = mat->GetNumberOfElements();
      G4double sum = 0.0;
      G4double norm= 0.0;
      for(G4int j=0; j<nelm; ++j) {
        G4int Z = (*elmVector)[j]->GetZasInt();
        G4double mass = fNistManager->GetAtomicMassAmu(Z)*CLHEP::amu_c2;
        G4int Z2 = Z*Z;
        sum += mass*Z2;
        norm += Z2;
      }
      effMass[i] = sum/norm;
    }
  }
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4WentzelVIRelModel::WentzelVIRelModelMutex);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
