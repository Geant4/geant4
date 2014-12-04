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
/*
 * G4ITManagerTest.cc
 *
 *  Created on: 5 mars 2014
 *      Author: kara
 */

#define MOLECULE_COUNTER_TESTING
#include "G4ITManager2.hh"
#include "G4Molecule.hh"
#include "G4OH.hh"
#include "G4Electron_aq.hh"

int main()
{
	G4ITManager2::Instance();
	G4Molecule* e_aq = new G4Molecule(G4Electron_aq::Definition());
        G4Molecule* OH = new G4Molecule(G4OH::Definition());

	// G4ITManager2::Instance()->PushMain(e_aq->BuildTrack(0,G4ThreeVector()));
    // G4ITManager2::Instance()->PushMain(OH->BuildTrack(0,G4ThreeVector()));

    G4ITManager2::Instance()->Push(e_aq->BuildTrack(0,G4ThreeVector()), ITMain);
    G4ITManager2::Instance()->Push(OH->BuildTrack(0,G4ThreeVector()), ITDelayed);

	ListGroup* listGroup = G4ITManager2::Instance()->GetListGroup();

	ListGroup::all_iterator it = listGroup->begin(ITMain);

	for( ; it != listGroup->end(ITMain) ; it ++)
	{
		G4cout << "in loop"  << G4endl;
		G4cout << GetIT(*it)->GetName() << G4endl;
	}


//        G4ITManager2::Instance()->Clear();

// Pop main deletes the tracks
	G4ITManager2::Instance()->PopMain(e_aq->GetTrack());
    G4ITManager2::Instance()->PopMain(OH->GetTrack());

	G4ITManager2::Instance()->Clear();

	G4cout << "test OK" << G4endl;
	return 0;
}



