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

#include "Test29VSecondaryGenerator.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicInteraction.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"


Test29VSecondaryGenerator::Test29VSecondaryGenerator(G4HadronicInteraction* hadi,
                                                     G4Material* mat):
  hInteraction(hadi), material(mat)
{
  G4cout<<"Test29VSecondaryGenerator::Constructor: material="<<material->GetName()<<G4endl;
  generatorName = "";
  const G4ElementVector* elv = material->GetElementVector();
  G4int Z = (G4int)(((*elv)[0])->GetZ() + 0.5);
  G4int N = (G4int)(((*elv)[0])->GetN() + 0.5);
  G4cout<<"Test29VSecondaryGenerator::Constructor: Nucleus with N="<<N<<", Z="<<Z<<G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  G4cout<<"Test29VSecondaryGenerator::Constr: Mass(targetNucleus)(MeV)="<<mass/MeV<<G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  G4cout << "Test29VSecondaryGenerator::Construct: Mass(IonTable)(MeV)="<<mass/MeV<<G4endl;
}

Test29VSecondaryGenerator::~Test29VSecondaryGenerator()
{}

// This member function starts the Generator from the HadronicModelCollection
G4HadFinalState* Test29VSecondaryGenerator::Secondaries(const G4Track& track)
{
  G4HadProjectile thePro(track); // One more HadronicPackageClass for the Projectile (?)
  if (hInteraction)
  {
    G4HadFinalState *result = hInteraction->ApplyYourself(thePro, targetNucleus);
    result->SetTrafoToLab(thePro.GetTrafoToLab()); // Some of Generators are not in LabSys?
    return result;
    //return hInteraction->ApplyYourself(track, targetNucleus); // Old fassion
  }
  else
  {
    G4cout<<"**Test29VSecondaryGenerator::Secondaries:Not initialized Interaction"<<G4endl;
    return 0; // The interaction is not initialized -> wrong initialization of the class
  }
}
