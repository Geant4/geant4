// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testMaterialPropertiesTable.cc,v 1.2 1999-12-15 14:50:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------

#include "G4ios.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"

void LoopUntilPressEnter();

int main()
{
	// Some Test Data

        const G4int NUMENTRIES = 32;

        G4double PPCKOV[NUMENTRIES] =
                  { 2.038E-9, 2.072E-9, 2.107E-9, 2.143E-9, 2.181E-9,
                    2.220E-9, 2.260E-9, 2.302E-9, 2.346E-9, 2.391E-9,
                    2.438E-9, 2.486E-9, 2.537E-9, 2.590E-9, 2.645E-9,
                    2.702E-9, 2.763E-9, 2.825E-9, 2.891E-9, 2.960E-9,
                    3.032E-9, 3.108E-9, 3.188E-9, 3.271E-9, 3.360E-9,
                    3.453E-9, 3.552E-9, 3.656E-9, 3.767E-9, 3.884E-9,
                    4.010E-9, 4.144E-9 };

        G4double EFFIC[NUMENTRIES] =
                 {  0.005,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
                    0.08, 0.09, 0.10, 0.115,0.13, 0.15, 0.16, 0.18,
                    0.195,0.22, 0.23, 0.24, 0.25, 0.255,0.26, 0.265,
                    0.26, 0.25, 0.24, 0.215,0.175,0.14, 0.085, 0.0 };


        G4double RINDEX[NUMENTRIES] =
                 {  1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33,
                    1.33, 1.33, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
                    1.34, 1.34, 1.35, 1.35, 1.35, 1.35, 1.35,
                    1.35, 1.35, 1.35, 1.35 };


        G4double ABSCO[NUMENTRIES] =
                 {  344.8,  408.2,  632.9,  917.4, 1234.6, 1388.9,
                    1515.2, 1724.1, 1886.8, 2000.0, 2631.6, 3571.4,
                    4545.5, 4761.9, 5263.2, 5263.2, 5555.6, 5263.2,
                    5263.2, 4761.9, 4545.5, 4166.7, 3703.7, 3333.3,
                    3000.0, 2850.0, 2700.0, 2450.0, 2200.0, 1950.0,
                    1750.0, 1450.0 };

	G4MaterialPropertiesTable myMPT;

//////////////////////////////////////////////////////////////////////////////
// Test storage and Retrieval
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\n\n\nTest AddProperty() and GetProperty()\n";
	G4cout << "------------------------------------\n\n"; 
	G4cout << "Store and Retrieve EFFIC and ABSCO\n\n";
        LoopUntilPressEnter();

	myMPT.AddProperty("EFFIC", PPCKOV, EFFIC, NUMENTRIES);

	G4MaterialPropertyVector *mpv;
	mpv = myMPT.GetProperty("EFFIC"); 

	G4cout << "EFFIC" << G4endl;
	G4cout << "-----" << G4endl;

	mpv->DumpVector();
	LoopUntilPressEnter();	

	myMPT.AddProperty("ABSCO", PPCKOV, ABSCO, NUMENTRIES);
	mpv = myMPT.GetProperty("ABSCO");

        G4cout << "ABSCO" << G4endl;
        G4cout << "-----" << G4endl;
        mpv->DumpVector();

	LoopUntilPressEnter();	

//////////////////////////////////////////////////////////////////////////////
// Test Generation of RINDEX entry
//////////////////////////////////////////////////////////////////////////////

	G4cout << "AddProperty RINDEX and Retrieve it again \n\n";
        LoopUntilPressEnter();

	myMPT.AddProperty("RINDEX", PPCKOV, RINDEX, NUMENTRIES);
	mpv = myMPT.GetProperty("RINDEX"); 

	G4cout << "RINDEX" << G4endl;
	G4cout << "------" << G4endl;

        mpv->DumpVector();

	LoopUntilPressEnter();	

	// Test Remove Property
	// --------------------
	G4cout << "\nTest RemoveProperty() -- Remove EFFIC\n";
	G4cout << "-------------------------------------\n\n";
	G4cout << "Dump table contents to observe absence of EFFIC\n\n";
        LoopUntilPressEnter();

	myMPT.RemoveProperty("EFFIC");	
	myMPT.DumpTable();
	
//////////////////////////////////////////////////////////////////////////////
// Test AddEntry 
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest AddEntry()\n";
	G4cout << "---------------\n\n";
	G4cout << "Add an element (3.166e-09, 1.34) to the Refraction Index \n";
	G4cout << "Property vector, check to see that it has been \n";
	G4cout << "inserted in its proper place. \n\n"; 
        LoopUntilPressEnter();

	mpv = myMPT.GetProperty("RINDEX");
	myMPT.AddEntry("RINDEX", 3.166e-09, 1.34);

        G4cout << "RINDEX" << G4endl;
        G4cout << "------" << G4endl;
	mpv->DumpVector();

	G4cout << "TESTING GetPhotonMomentum()" << G4endl;
	G4cout << "---------------------------" << G4endl;

	for (int i=0; i < NUMENTRIES; i++) { 
		G4cerr << "Photon Momentum for " << RINDEX[i] << " "
		     << mpv->GetPhotonMomentum(RINDEX[i])	
		     << G4endl;
	}	

        LoopUntilPressEnter();

//////////////////////////////////////////////////////////////////////////////
// Test RemoveEntry 
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest RemoveEntry()\n";
	G4cout << "------------------\n\n";
	G4cout << "Remove an element (3.166e-09, 1.34) from the Refraction \n";
	G4cout << "Index Property vector, check to see that it has been \n";
	G4cout << "removed from its proper place \n\n";
        LoopUntilPressEnter();

        mpv = myMPT.GetProperty("RINDEX");
        myMPT.RemoveEntry("RINDEX", 3.166e-09);

        G4cout << "RINDEX" << G4endl;
        G4cout << "------" << G4endl;
        mpv->DumpVector();
        LoopUntilPressEnter();

        myMPT.AddProperty("EFFIC", PPCKOV, EFFIC, NUMENTRIES);

//////////////////////////////////////////////////////////////////////////////
// Test Assignment operator for MaterialPropertyVector  
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest Assignment Operator for MaterialPropertyVector\n";
	G4cout << "--------------------------------------------------\n\n";

        LoopUntilPressEnter();

	G4MaterialPropertyVector testMPV1;
	G4MaterialPropertyVector testMPV2;
	
	testMPV1 = *myMPT.GetProperty("EFFIC");

	G4cout << "Vector1\n\n";

	testMPV1.DumpVector();
        LoopUntilPressEnter();

	testMPV2 = *myMPT.GetProperty("RINDEX");

	G4cout << "Vector2\n\n";

	testMPV2.DumpVector();
        LoopUntilPressEnter();

	testMPV1 = testMPV2;

        G4cout << "Vector1\n\n";

        testMPV1.DumpVector();
        LoopUntilPressEnter();

//////////////////////////////////////////////////////////////////////////////
// Test Copy Constructor for MaterialPropertyVector  
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest Copy Constructor for MaterialPropertyVector\n";
	G4cout << "-----------------------------------------------\n\n";
	G4cout << "Copying contents of Vector2 into Vector3\n\n";
        LoopUntilPressEnter();

	G4MaterialPropertyVector testMPV3(testMPV2);
        G4cout << "Vector3\n\n";

        testMPV2.DumpVector();

        G4cout << "MinPM " << testMPV2.GetMinPhotonMomentum() << G4endl;
	G4cout << "MaxPM " << testMPV2.GetMaxPhotonMomentum() << G4endl;
        G4cout << "\n\n";
        LoopUntilPressEnter();

//////////////////////////////////////////////////////////////////////////////
// Test Assignment operator for MaterialPropertiesTable 
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest Assignment operator for MaterialPropertiesTable\n";
	G4cout << "---------------------------------------------------\n\n";
	G4cout << "Assign Table2 the contents of the current"; 
	G4cout << "MaterialPropertiesTable\n\n";
        LoopUntilPressEnter();

	G4MaterialPropertiesTable myMPT2; 
	myMPT2 = myMPT;

	G4cout << "Table2\n\n";

	myMPT2.DumpTable();

//////////////////////////////////////////////////////////////////////////////
// Test Copy Constructor for MaterialPropertiesTable 
//////////////////////////////////////////////////////////////////////////////

	G4cout << "\nTest Copy Constructor for MaterialPropertiesTable\n";
	G4cout << "------------------------------------------------\n\n";
	G4cout << "Assign Table3 the contents of Table2\n\n";
        LoopUntilPressEnter();

	G4MaterialPropertiesTable myMPT3(myMPT2);

	G4cout << "Table3\n\n";

	myMPT3.DumpTable();

	G4cout << "\n\n\n<END OF TEST>\n\n\n";

	return EXIT_SUCCESS;
}

// LoopUntilPressEnter
// -------------------
//
void LoopUntilPressEnter()
{
        char ch;
	G4cout << "Press <Enter> to continue ... ";
        while ( G4cin.get(ch) )
        {
                if (ch == '\n') break;
        }       
	G4cout << G4endl;
}

