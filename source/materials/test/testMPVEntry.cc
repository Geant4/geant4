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
//
// $Id: testMPVEntry.cc,v 1.4 2001-07-17 15:54:44 verderi Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "G4MPVEntry.hh"

G4MPVEntry TestStorage (void);
void TestLessThan(G4MPVEntry& e1, G4MPVEntry& e2);
void TestLogicalEquality(G4MPVEntry& e1, G4MPVEntry& e2);	

int main()
{

	G4MPVEntry Entry1(TestStorage());
	G4MPVEntry Entry2(TestStorage());
	TestLessThan(Entry1, Entry2);
	TestLogicalEquality(Entry1, Entry2);	

	return EXIT_SUCCESS;
}


G4MPVEntry TestStorage (void)
{
	G4double pm, pv;
	
 	G4cout << "Input photon momentum for G4MPVEntry:  " << G4endl;
        G4cin >> pm;
        G4cout << "Input a property value for this photon momentum: " << G4endl;
        G4cin >> pv;
        
        G4MPVEntry newEntry(pm, pv);

        G4cout << "This is your entry: " << G4endl;
        newEntry.DumpEntry();
	G4cout << G4endl;
	G4cout << G4endl;

	return newEntry;
}

void TestLessThan(G4MPVEntry& e1, G4MPVEntry& e2)
{
	if (e1 < e2)
		G4cout << "first entry less than second." << G4endl;
	else { 
		G4cout << "first entry greater than or equal to second" << G4endl;
	}
}

void TestLogicalEquality(G4MPVEntry& e1, G4MPVEntry& e2)
{
	if (e1 == e2)
		G4cout << "the two entries are equal." << G4endl;
	else  
		G4cout << "the two entries are NOT equal." << G4endl;

}




