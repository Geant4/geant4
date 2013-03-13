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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RadioactiveDecayPhysics.hh"

#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4RadioactiveDecayPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics::G4RadioactiveDecayPhysics(G4int)
  :  G4VPhysicsConstructor("G4RadioactiveDecay"), theRadioactiveDecay(0)
{}

G4RadioactiveDecayPhysics::G4RadioactiveDecayPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name), theRadioactiveDecay(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics::~G4RadioactiveDecayPhysics()
{
  delete theRadioactiveDecay;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics::ConstructProcess()
{
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(new G4RadioactiveDecay(), G4GenericIon::GenericIon());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

