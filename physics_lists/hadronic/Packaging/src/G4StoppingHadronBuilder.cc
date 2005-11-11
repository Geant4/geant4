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
// $Id: G4StoppingHadronBuilder.cc,v 1.2 2005-11-11 11:13:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4StoppingHadronBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard and add mu-
//
//----------------------------------------------------------------------------
//

#include "G4StoppingHadronBuilder.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"

G4StoppingHadronBuilder::G4StoppingHadronBuilder(): wasActivated(false)
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
