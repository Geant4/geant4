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

#include "BiasedRDPhysics.hh"

#include "G4Radioactivation.hh"
#include "G4GenericIon.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"
#include "G4EmParameters.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(BiasedRDPhysics);


BiasedRDPhysics::BiasedRDPhysics(G4int)
 : G4VPhysicsConstructor("G4Radioactivation")
{
  G4EmParameters::Instance()->AddPhysics("World","G4Radioactivation");
  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetStoreICLevelData(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                       /std::log(2.));
}

BiasedRDPhysics::BiasedRDPhysics(const G4String&)
 : BiasedRDPhysics(0)
{}


BiasedRDPhysics::~BiasedRDPhysics()
{}


void BiasedRDPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
}


void BiasedRDPhysics::ConstructProcess()
{
  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if(!ad) {
    G4EmParameters::Instance()->SetAugerCascade(true);
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    ad->InitialiseAtomicDeexcitation();
  }

  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(new G4Radioactivation(), G4GenericIon::GenericIon());
}

