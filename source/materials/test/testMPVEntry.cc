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
// $Id: testMPVEntry.cc,v 1.5 2006-06-29 19:13:31 gunter Exp $
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




