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
  material(mat)
{
  G4cout << "New generator and material= " << material->GetName() << G4endl;
  generatorName = "";
  const G4ElementVector* ev = material->GetElementVector();
  G4int Z = (G4int)(((*ev)[0])->GetZ() + 0.5);
  G4int N = (G4int)(((*ev)[0])->GetN() + 0.5);
  G4cout << "Nucleus with N= " << N << "  Z= " << Z << G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  G4cout << "Mass from targetNucleus(MeV)= " << mass/MeV << G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  G4cout << "Mass from IonTable(MeV)=      " << mass/MeV << G4endl;	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30VSecondaryGenerator::~Test30VSecondaryGenerator()
{
  //  if(hInteraction) delete hInteraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* Test30VSecondaryGenerator::Secondaries(const G4Track& track)
{
 if (hInteraction) {
  
   return hInteraction->ApplyYourself(track, targetNucleus);

 } else {
   theParticleChange.Initialize(track);
   return &theParticleChange;
 }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






