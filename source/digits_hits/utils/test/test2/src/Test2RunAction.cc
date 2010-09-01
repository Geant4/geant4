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
// $Id: Test2RunAction.cc,v 1.2 2010-09-01 08:03:10 akimura Exp $
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

  const Test2Run* theRun = (const Test2Run*)aRun;
  
  G4cout << "############################################################" << G4endl;
  G4cout << " Run Summary - Number of events : " << theRun->GetNumberOfEvent() << G4endl;
  G4cout << "############################################################" << G4endl;

  // output of primitive scorers
  G4cout << "                   == Primitive Scorers ==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRun->GetTotal(0),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRun->GetTotal(1),"Length")
	 << "   nStep " << theRun->GetTotal(2) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRun->GetTotal(3),"Length")
	 << "   nStep " << theRun->GetTotal(4) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRun->GetTotal(5),"Length")
	 << "   nStep " << theRun->GetTotal(6) << G4endl;
  G4cout << G4endl;

  // output of traditional sensitive manager
  G4cout << "           == Traditional Sensitive Detector ==" << G4endl;
  G4cout << "Total energy deposition in phantom : "
	 << G4BestUnit(theRun->GetSDTotalEdep(),"Energy") << G4endl;
  G4cout << "Gamma -       track length " << G4BestUnit(theRun->GetSDTrackLengthGamma(),"Length") << G4endl;
  //	 << "   nStep " << theRun->GetTotal(2) << G4endl;
  G4cout << "Electron -    track length " << G4BestUnit(theRun->GetSDTrackLengthElec(),"Length") << G4endl;
    //<< "   nStep " << theRun->GetTotal(4) << G4endl;
  G4cout << "Positron -    track length " << G4BestUnit(theRun->GetSDTrackLengthPosi(),"Length") << G4endl;
    //	 << "   nStep " << theRun->GetTotal(6) << G4endl;
  G4cout << G4endl;
  
}

