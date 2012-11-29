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
// ClassName:   G4LHEPStoppingPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard and add mu-
// 16.10.2012 A.Ribon: deprecated class, to be removed in G4 version 10.
//
//----------------------------------------------------------------------------
//

#include "G4LHEPStoppingPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MuonMinus.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4MuonMinusCaptureAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4HadronicDeprecate.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4LHEPStoppingPhysics);


G4LHEPStoppingPhysics::G4LHEPStoppingPhysics(G4int ver)
  : G4VPhysicsConstructor("LHEP Stopping")
 , muProcess(0), piProcess(0), kProcess(0), apProcess(0), anProcess(0)
    , verbose(ver), wasActivated(false)
{
  G4HadronicDeprecate("G4LHEPStoppingPhysics");
}

G4LHEPStoppingPhysics::G4LHEPStoppingPhysics(const G4String& nam, G4int ver)
: G4VPhysicsConstructor(nam), muProcess(0), piProcess(0), kProcess(0), apProcess(0), anProcess(0)
, verbose(ver), wasActivated(false)
{
  G4HadronicDeprecate("G4LHEPStoppingPhysics");
}

G4LHEPStoppingPhysics::~G4LHEPStoppingPhysics()
{
  if(wasActivated) {
    delete muProcess;
    delete piProcess;
    delete kProcess;
    delete apProcess;
    delete anProcess;
  }
}

void G4LHEPStoppingPhysics::ConstructParticle()
{
// G4cout << "G4QStoppingPhysics::ConstructParticle" << G4endl;
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

}

void G4LHEPStoppingPhysics::ConstructProcess()
{
  if(wasActivated) return;
  G4ProcessManager * aProcMan = 0;
  wasActivated=true;
  //G4cout << " adding stopping hadron Physics" << G4endl;

  // Muon Minus Physics
  aProcMan = G4MuonMinus::MuonMinus()->GetProcessManager();
  muProcess = new G4MuonMinusCaptureAtRest();
  aProcMan->AddRestProcess(muProcess);

  // PionMinus
  aProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  piProcess = new G4PionMinusAbsorptionAtRest();
  aProcMan->AddRestProcess(piProcess);

  // KaonMinus
  aProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  kProcess = new G4KaonMinusAbsorption();
  aProcMan->AddRestProcess(kProcess);

  // anti-Proton
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  apProcess = new G4AntiProtonAnnihilationAtRest();
  aProcMan->AddRestProcess(apProcess);

  // AntiNeutron
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  anProcess = new G4AntiNeutronAnnihilationAtRest();
  aProcMan->AddRestProcess(anProcess);

}
