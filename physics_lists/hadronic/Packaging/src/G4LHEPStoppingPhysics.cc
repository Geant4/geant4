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
// $Id: G4LHEPStoppingPhysics.cc,v 1.1 2006-05-04 16:48:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LHEPStoppingPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard and add mu-
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

G4LHEPStoppingPhysics::G4LHEPStoppingPhysics(const G4String& nam, G4int ver)
  : G4VPhysicsConstructor(nam), verbose(ver), wasActivated(false)
{}

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
