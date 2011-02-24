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
// $Id: Test2RunAction.cc,v 1.4 2010-11-03 08:48:57 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Test2RunAction.hh"
#include "Test2Run.hh"

#include "G4ios.hh"
#include "G4UnitsTable.hh"

Test2RunAction::Test2RunAction()
{;}

Test2RunAction::~Test2RunAction()
{;}

G4Run* Test2RunAction::GenerateRun() {
  return new Test2Run; 
}

void Test2RunAction::BeginOfRunAction(const G4Run*)
{;}

void Test2RunAction::EndOfRunAction(const G4Run* aRun) {

  Test2Run* theRun = (Test2Run*)aRun;

  theRun->DumpQuantitiesToFile();
  
  G4cout << "############################################################" << G4endl;
  G4cout << " Run Summary - Number of events : " << theRun->GetNumberOfEvent() << G4endl;
  G4cout << "############################################################" << G4endl;


  // output of traditional sensitive manager
  const Test2RunSD* theRunSD = theRun->GetMassRunSD();
  G4cout << "           == Traditional Sensitive Detector in MassWorld==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRunSD->GetTotal(0),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRunSD->GetTotal(1),"Length")
  	 << "   nStep " << theRunSD->GetTotal(4) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRunSD->GetTotal(2),"Length")
	 << "   nStep " << theRunSD->GetTotal(5) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRunSD->GetTotal(3),"Length")
    	 << "   nStep " << theRunSD->GetTotal(6) << G4endl;
  G4cout << G4endl;

  // output of primitive scorers
  const Test2RunPS* theRunPS = theRun->GetMassRunPS();
  G4cout << "                   == Primitive Scorers in MassWorld==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRunPS->GetTotal(0),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRunPS->GetTotal(1),"Length")
	 << "   nStep " << theRunPS->GetTotal(4) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRunPS->GetTotal(2),"Length")
	 << "   nStep " << theRunPS->GetTotal(5) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRunPS->GetTotal(3),"Length")
	 << "   nStep " << theRunPS->GetTotal(6) << G4endl;
  G4cout << G4endl;

  // output of traditional sensitive manager
  const Test2RunSD* theRunSDp = theRun->GetParaRunSD();
  G4cout << "           == Traditional Sensitive Detector in ParallelWorld==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRunSDp->GetTotal(0),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRunSDp->GetTotal(1),"Length")
  	 << "   nStep " << theRunSDp->GetTotal(4) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRunSDp->GetTotal(2),"Length")
	 << "   nStep " << theRunSDp->GetTotal(5) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRunSDp->GetTotal(3),"Length")
    	 << "   nStep " << theRunSDp->GetTotal(6) << G4endl;
  G4cout << G4endl;

  // output of primitive scorers
  theRunPS = theRun->GetParaRunPS();
  G4cout << "                   == Primitive Scorers in ParallelWorld==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRunPS->GetTotal(0),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRunPS->GetTotal(1),"Length")
	 << "   nStep " << theRunPS->GetTotal(4) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRunPS->GetTotal(2),"Length")
	 << "   nStep " << theRunPS->GetTotal(5) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRunPS->GetTotal(3),"Length")
	 << "   nStep " << theRunPS->GetTotal(6) << G4endl;
  G4cout << G4endl;
  
}

