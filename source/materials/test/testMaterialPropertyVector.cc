// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testMaterialPropertyVector.cc,v 1.1 1999-01-08 16:32:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------

#include "G4MaterialPropertyVector.hh"

void LoopUntilPressEnter();

int main ()
{
        const G4int NUMENTRIES = 32;
	G4double aPhotonMomentum, anAbsorptionCoefficient;

        G4double PPCKOV[NUMENTRIES] =
                  { 2.038E-9, 2.072E-9, 2.107E-9, 2.143E-9, 2.181E-9,
                    2.220E-9, 2.260E-9, 2.302E-9, 2.346E-9, 2.391E-9,
                    2.438E-9, 2.486E-9, 2.537E-9, 2.590E-9, 2.645E-9,
                    2.702E-9, 2.763E-9, 2.825E-9, 2.891E-9, 2.960E-9,
                    3.032E-9, 3.108E-9, 3.188E-9, 3.271E-9, 3.360E-9,
                    3.453E-9, 3.552E-9, 3.656E-9, 3.767E-9, 3.884E-9,
                    4.010E-9, 4.144E-9 };

        G4double ABSCO[NUMENTRIES] =
                 {  344.8,  408.2,  632.9,  917.4, 1234.6, 1388.9,
                    1515.2, 1724.1, 1886.8, 2000.0, 2631.6, 3571.4,
                    4545.5, 4761.9, 5263.2, 5263.2, 5555.6, 5263.2,
                    5263.2, 4761.9, 4545.5, 4166.7, 3703.7, 3333.3,
                    3000.0, 2850.0, 2700.0, 2450.0, 2200.0, 1950.0,
                    1750.0, 1450.0 };


	// Test Vector creation
	// --------------------
	G4cout << "Test Vector creation" << endl;
	G4cout << "--------------------" << endl << endl;

	G4MaterialPropertyVector absco(PPCKOV, ABSCO, NUMENTRIES);
	absco.DumpVector();
	LoopUntilPressEnter();
		
	// Test GetProperty
	// ----------------	
	G4cout << "Test GetProperty" << endl;
	G4cout << "----------------" << endl << endl; 
        G4cout << "Input a photon momentum value for which you wish" << endl;
        G4cout << "to find the absorption coefficient:  ";

        cin >> aPhotonMomentum;
	G4cout << endl;

	anAbsorptionCoefficient = absco.GetProperty(aPhotonMomentum);

	G4cout << "The absorption coefficient is:  " 
	     << anAbsorptionCoefficient  
	     << endl;
        LoopUntilPressEnter();

	// Test AddElement
	// ---------------
	G4cout << "Test AddElement" << endl;
	G4cout << "---------------" << endl << endl;
	G4cout << "Add the entry just created" << endl;
	LoopUntilPressEnter();

	absco.AddElement(aPhotonMomentum, anAbsorptionCoefficient); 
	absco.DumpVector();

	// Test RemoveElement
	// ------------------
	G4cout << "Test RemoveElement" << endl;
        G4cout << "------------------" << endl << endl;
        G4cout << "Input a photon momentum value for which you wish" << endl;
        G4cout << "to remove the OPVEntry:  ";

        cin >> aPhotonMomentum;
	G4cout << endl;

	absco.RemoveElement(aPhotonMomentum);
	absco.DumpVector(); 

	// Test the iterator
	// -----------------
	G4cout << "Test the iterator" << endl;
	G4cout << "-----------------" << endl << endl;
        LoopUntilPressEnter();

	absco.ResetIterator();

	while (++absco)
	{
		G4cout << absco.GetPhotonMomentum() << "\t"
		     << absco.GetProperty() << endl;
	}
	G4cout << "\n\n <END OF TEST>" << endl;

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
