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



