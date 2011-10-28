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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test30VSecondaryGenerator -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test30VSecondaryGenerator.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicInteraction.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"
#include "G4NucleiProperties.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30VSecondaryGenerator::Test30VSecondaryGenerator(G4HadronicInteraction* hadi,
  G4Material* mat):
  hInteraction(hadi),
  material(mat),
  targetNucleus(mat),
  targetN(0.0)
{
  elm = material->GetElement(0);
  targetZ = elm->GetZ();
  G4cout << "New generator for material " << material->GetName() 
	 << " Nelm= " <<  material->GetNumberOfElements() 
	 << " Nmat= " <<  material->GetNumberOfMaterials() 
	 << " Target element: " << elm->GetName() 
	 << G4endl;
  generatorName = "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30VSecondaryGenerator::~Test30VSecondaryGenerator()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test30VSecondaryGenerator::SetA(G4int A) 
{
  targetN = G4double(A);
  G4int Z = (G4int)(targetZ + 0.5);
  G4int N = (G4int)(elm->GetN() + 0.5);
  if(targetN > 0.0) N = (G4int)(targetN + 0.5);
  G4cout << "Nucleus with N= " << N << "  Z= " << Z << "  A(amu)= " << elm->GetN();
  if(targetN == 0.0) G4cout << "  Natural abandances";
  G4cout << G4endl;
  G4double mass = G4NucleiProperties::GetNuclearMass(G4double(N), targetZ);
  G4cout << "Mass from targetNucleus(GeV)= " << mass/GeV << G4endl;
  G4double mass1 = 
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  G4cout << "Mass from IonTable(GeV)=      " << mass1/GeV << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double Test30VSecondaryGenerator::GetMass()
{
  G4double currentN = targetN;
  if(targetN == 0.0) {
    currentN = elm->GetN();
    G4IsotopeVector* isoVector = elm->GetIsotopeVector();
    G4int nIsoPerElement = elm->GetNumberOfIsotopes();
    currentN = G4double((*isoVector)[0]->GetN());
    if(nIsoPerElement > 1) {
      G4double* abundVector = elm->GetRelativeAbundanceVector();
      G4double y = G4UniformRand();
      for(G4int j=0; j<nIsoPerElement; j++) {
	y -= abundVector[j];
	if(y <= 0.0) {
	  currentN = G4double((*isoVector)[j]->GetN());
	  break;
	}
      }
    }
  }
  targetNucleus.SetParameters(currentN, targetZ);
  G4double mass = G4NucleiProperties::GetNuclearMass(currentN, targetZ);
  return mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* Test30VSecondaryGenerator::Secondaries(const G4Track& track)
{
  G4HadProjectile thePro(track);
  if (hInteraction) {

    G4HadFinalState *result = hInteraction->ApplyYourself(thePro, targetNucleus);
    result->SetTrafoToLab(thePro.GetTrafoToLab());
    return result;

 } else {
    return 0;
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
