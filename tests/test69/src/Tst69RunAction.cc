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
// $Id$
//

#include "Tst69RunAction.hh"
#include "Tst69INCLXXTallyAnalysis.hh"
#ifdef TEST69_HAS_ROOT
#include "Tst69INCLXXTallyROOT.hh"
#endif

#include "G4Run.hh"

#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"

#include "G4INCLXXInterface.hh"
#include "G4INCLXXInterfaceStore.hh"

#include "G4AblaInterface.hh"
#include "Randomize.hh"

#include <cstdlib>
#include <vector>

Tst69RunAction::Tst69RunAction(const char * const physList) {
  const char *tallyOption = getenv("TEST69_TALLY");
  G4String g4TallyOption;
  if(tallyOption)
    g4TallyOption = tallyOption;
  if(g4TallyOption=="analysis") {
    theTally = new Tst69INCLXXTallyAnalysis(physList);
    G4cout << "TallyAnalysis selected" << G4endl;
  }
#ifdef TEST69_HAS_ROOT
  else if(g4TallyOption=="root") {
    const int tID = G4Threading::G4GetThreadId();
    if(tID==-2) {
      theTally = new Tst69INCLXXTallyROOT(physList);
      G4cout << "TallyROOT selected" << G4endl;
    } else if(tID==0) {
      theTally = new Tst69INCLXXTallyROOT(physList);
      G4cout << "TallyROOT selected on worker thread " << tID << G4endl;
    } else {
      theTally = NULL;
      G4cout << "WARNING: TallyROOT requested in MT mode, but ROOT is not thread-safe." << G4endl;
      if(tID==-1)
        G4cout << "No tally selected on master thread" << G4endl;
      else
        G4cout << "No tally selected on worker thread " << tID << G4endl;
    }
  }
#endif
  else {
    theTally = NULL;
    G4cout << "no tally selected" << G4endl;
  }
}

Tst69RunAction::~Tst69RunAction()
{
  delete theTally;
}

void Tst69RunAction::BeginOfRunAction(const G4Run* )
{
  if(getenv("TEST69_USE_ABLA")) {
    std::vector<G4HadronicInteraction *> interactions = G4HadronicInteractionRegistry::Instance()
      ->FindAllModels(G4INCLXXInterfaceStore::GetInstance()->getINCLXXVersionName());
    for(std::vector<G4HadronicInteraction *>::const_iterator iInter=interactions.begin(), e=interactions.end();
        iInter!=e; ++iInter) {
      G4INCLXXInterface *theINCLInterface = static_cast<G4INCLXXInterface*>(*iInter);
      if(theINCLInterface) {
        G4HadronicInteraction *interaction = G4HadronicInteractionRegistry::Instance()->FindModel("ABLA");
        G4AblaInterface *theAblaInterface = static_cast<G4AblaInterface*>(interaction);
        if(!theAblaInterface)
          theAblaInterface = new G4AblaInterface;
        G4cout << "Coupling INCLXX to ABLA" << G4endl;
        theINCLInterface->SetDeExcitation(theAblaInterface);
      }
    }
  }

  // set the INCL++ tally object if necessary
  if(theTally) {
    G4INCLXXInterfaceStore *theStore = G4INCLXXInterfaceStore::GetInstance();
    if(theStore) {
      G4cout << "Activating tally for INCLXX" << G4endl;
      theStore->SetTally(theTally);
      theTally->Open();
    }
  }
}

void Tst69RunAction::EndOfRunAction(const G4Run*) {
  if(theTally)
    theTally->Close();
}

