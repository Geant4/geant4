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

#include "Test49VSecondaryGenerator.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicInteraction.hh"
#include "G4ElementVector.hh"
#include "G4Material.hh"


Test49VSecondaryGenerator::Test49VSecondaryGenerator(G4HadronicInteraction* hadi,
                                                     G4Material* mat):
  hInteraction(hadi), material(mat)
{
  G4cout<<"Test49VSecondaryGenerator::Constructor: material="<<material->GetName()<<G4endl;
  generatorName = "";
  const G4ElementVector* elv = material->GetElementVector();
  G4int Z = (G4int)(((*elv)[0])->GetZ() + 0.5);
  G4int N = (G4int)(((*elv)[0])->GetN() + 0.5);
  G4cout<<"Test49VSecondaryGenerator::Constructor: Nucleus with N="<<N<<", Z="<<Z<<G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  G4cout<<"Test49VSecondaryGenerator::Constr: Mass(targetNucleus)(MeV)="<<mass/MeV<<G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  G4cout << "Test49VSecondaryGenerator::Construct: Mass(IonTable)(MeV)="<<mass/MeV<<G4endl;
}

Test49VSecondaryGenerator::~Test49VSecondaryGenerator()
{}

// This member function starts the Generator from the HadronicModelCollection
G4HadFinalState* Test49VSecondaryGenerator::Secondaries(const G4Track& track)
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
    G4cout<<"**Test49VSecondaryGenerator::Secondaries:Not initialized Interaction"<<G4endl;
    return 0; // The interaction is not initialized -> wrong initialization of the class
  }
}
