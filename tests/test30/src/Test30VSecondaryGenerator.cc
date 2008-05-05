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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30VSecondaryGenerator::Test30VSecondaryGenerator(G4HadronicInteraction* hadi,
  G4Material* mat):
  hInteraction(hadi),
  material(mat),
  targetN(0)
{
  G4cout << "New generator and material= " << material->GetName() 
	 << " Nelm= " <<  material->GetNumberOfElements() 
	 << " Nmat= " <<  material->GetNumberOfMaterials() 
	 << G4endl;
  generatorName = "";
  const G4ElementVector* ev = material->GetElementVector();
  G4int Z = (G4int)(((*ev)[0])->GetZ() + 0.5);
  G4int N = (G4int)(((*ev)[0])->GetN() + 0.5);
  if(targetN > 0) N = targetN;
  G4cout << "Nucleus with N= " << N << "  Z= " 
	 << Z << "  " << (*ev)[0]->GetN()<< G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  G4cout << "Mass from targetNucleus(MeV)= " << mass/MeV << G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  G4cout << "Mass from IonTable(MeV)=      " << mass/MeV << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30VSecondaryGenerator::~Test30VSecondaryGenerator()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* Test30VSecondaryGenerator::Secondaries(const G4Track& track)
{
  G4HadProjectile thePro(track);
  if (hInteraction) {

    G4HadFinalState *result = hInteraction->ApplyYourself(thePro, targetNucleus);
    result->SetTrafoToLab(thePro.GetTrafoToLab());
    return result;
//   return hInteraction->ApplyYourself(track, targetNucleus);

 } else {
    return 0;
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






