// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testMPVEntry.cc,v 1.1 1999-01-08 16:32:10 gunter Exp $
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
	
 	G4cout << "Input photon momentum for G4MPVEntry:  ";
        cin >> pm;
        G4cout << "Input a property value for this photon momentum: ";
        cin >> pv;
        
        G4MPVEntry newEntry(pm, pv);

        G4cout << "This is your entry: ";
        newEntry.DumpEntry();
	G4cout << endl;
	G4cout << endl;

	return newEntry;
}

void TestLessThan(G4MPVEntry& e1, G4MPVEntry& e2)
{
	if (e1 < e2)
		G4cout << "first entry less than second." << endl;
	else { 
		G4cout << "first entry greater than or equal to second" << endl;
	}
}

void TestLogicalEquality(G4MPVEntry& e1, G4MPVEntry& e2)
{
	if (e1 == e2)
		G4cout << "the two entries are equal." << endl;
	else  
		G4cout << "the two entries are NOT equal." << endl;

}
