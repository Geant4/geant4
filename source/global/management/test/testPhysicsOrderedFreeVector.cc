// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testPhysicsOrderedFreeVector.cc,v 1.1 1999-01-07 16:09:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "g4templates.hh"
#include "G4PhysicsOrderedFreeVector.hh"

void LoopUntilPressEnter();

int main ()
{
        const G4int NUMENTRIES = 32;
	G4double anEnergy, aValue;

        G4double PPCKOV[NUMENTRIES] =
                  { 2.038E-9, 2.072E-9, 2.107E-9, 2.143E-9, 2.181E-9,
                    2.220E-9, 2.260E-9, 2.302E-9, 2.346E-9, 2.391E-9,
                    2.438E-9, 2.486E-9, 2.537E-9, 2.590E-9, 2.645E-9,
                    2.702E-9, 2.763E-9, 2.825E-9, 2.891E-9, 2.960E-9,
                    3.032E-9, 3.108E-9, 3.188E-9, 3.271E-9, 3.360E-9,
                    3.453E-9, 3.552E-9, 3.656E-9, 3.767E-9, 3.884E-9,
                    4.010E-9, 4.144E-9 };

        G4double RINDEX[NUMENTRIES] =
                 {  1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33,
                    1.33, 1.33, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.35, 1.35, 1.35, 1.35, 1.35,
                    1.35, 1.35, 1.35, 1.35 };

	// Test Vector creation
	// --------------------
	G4cout << "Test Vector creation" << endl;
	G4cout << "--------------------" << endl << endl; 

	G4PhysicsOrderedFreeVector aVector(PPCKOV, RINDEX, NUMENTRIES);
	aVector.DumpValues();
	LoopUntilPressEnter();

	// Test GetEnergy
	// --------------
	G4cout << "Test GetEnergy" << endl;
	G4cout << "--------------" << endl;
	G4cout << "Input a value within the vector range for which you" << endl;
	G4cout << "wish to find the corresponding energy:  " << endl;
	cin >> aValue;

	anEnergy = aVector.GetEnergy(aValue);
	G4cout << "The corresponding energy is " << anEnergy << endl;

	// Test GetMaxValue 
	// ----------------	
	G4cout << "Test GetMaxValue" << endl;
	G4cout << "----------------" << endl << endl; 
        LoopUntilPressEnter();

	aValue = aVector.GetMaxValue();

	G4cout << "The Max Value is:  " << aValue << endl;
	LoopUntilPressEnter();

	// Test GetMinValue 
	// ----------------
	G4cout << "Test GetMinValue" << endl;
	G4cout << "----------------" << endl << endl;
 
	aValue = aVector.GetMinValue();

	G4cout << "The Max Value is:  " << aValue << endl;
	LoopUntilPressEnter();

	// Test GetMaxLowEdgeEnergy 
	// ------------------------
	G4cout << "Test GetMaxLowEdgeEnergy" << endl;
	G4cout << "------------------------" << endl << endl;
 
	anEnergy = aVector.GetMaxLowEdgeEnergy();

	G4cout << "The Max Value is:  " << anEnergy << endl;
	LoopUntilPressEnter();

	// Test GetMinLowEdgeEnergy 
	// ------------------------
	G4cout << "Test GetMinLowEdgeEnergy" << endl;
	G4cout << "------------------------" << endl << endl;
 
	anEnergy = aVector.GetMinLowEdgeEnergy();

	G4cout << "The Max Value is:  " << anEnergy << endl;

        return EXIT_SUCCESS;
}

// LoopUntilPressEnter
// -------------------
//
void LoopUntilPressEnter()
{
        char ch;
        G4cout << "Press <Enter> to continue ... ";
        while ( cin.get(ch) )
        {
                if (ch == '\n') break;
        }
        G4cout << endl;
}
