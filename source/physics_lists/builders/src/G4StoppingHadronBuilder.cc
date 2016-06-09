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
// ClassName:   G4StoppingHadronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard and add mu-
//----------------------------------------------------------------------------
//

#include "G4StoppingHadronBuilder.hh"
#include "G4MuonMinusCaptureAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"

G4StoppingHadronBuilder::G4StoppingHadronBuilder():
 theMuonMinusAbsorption(0),thePionMinusAbsorption(0),
 theKaonMinusAbsorption(0), theAntiProtonAnnihilation(0),
 theAntiNeutronAnnihilation(0),
 wasActivated(false)
{}

G4StoppingHadronBuilder::~G4StoppingHadronBuilder()
{
  if(wasActivated) {
    delete theMuonMinusAbsorption;
    delete thePionMinusAbsorption;
    delete theKaonMinusAbsorption;
    delete theAntiProtonAnnihilation;
    delete theAntiNeutronAnnihilation;
  }
}

void G4StoppingHadronBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  wasActivated=true;
//G4cout << " adding stopping hadron Physics" << G4endl;

  // Muon Minus Physics
  aProcMan = G4MuonMinus::MuonMinus()->GetProcessManager();
  theMuonMinusAbsorption = new G4MuonMinusCaptureAtRest();
  aProcMan->AddRestProcess(theMuonMinusAbsorption);

  // PionMinus
  aProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  thePionMinusAbsorption = new G4PionMinusAbsorptionAtRest();
  aProcMan->AddRestProcess(thePionMinusAbsorption);

  // KaonMinus
  aProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  theKaonMinusAbsorption = new G4KaonMinusAbsorption();
  aProcMan->AddRestProcess(theKaonMinusAbsorption);

  // anti-Proton
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  theAntiProtonAnnihilation = new G4AntiProtonAnnihilationAtRest();
  aProcMan->AddRestProcess(theAntiProtonAnnihilation);

  // AntiNeutron
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  theAntiNeutronAnnihilation = new G4AntiNeutronAnnihilationAtRest();
  aProcMan->AddRestProcess(theAntiNeutronAnnihilation);

}
