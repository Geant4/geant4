// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testAdditionOf_MPT_toMaterial.cc,v 1.2 1999-11-11 15:36:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"

#include "G4ios.hh"

//void G4RWTPtrOrderedVector<G4MPVEntry>::clearAndDestroy(void);
//G4MPVEntry* G4RWTPtrOrderedVector<G4MPVEntry>::find(const G4MPVEntry*) const;
//void G4RWTPtrOrderedVector<G4DecayChannel>::clearAndDestroy(void);

G4RWTPtrOrderedVector<G4MPVEntry> lof;
#include "G4DecayChannel.hh"
G4RWTPtrOrderedVector<G4DecayChannel> lof3;
void lof2() {
  lof.clearAndDestroy();
  lof3.clearAndDestroy();
  lof.find(0);
}


void LoopUntilPressEnter();

int main()
{
        //--------- Material definition ---------

	G4String name;
	G4String symbol;
	G4int iz, n, nel;
	G4double z, a, density;

	a = 1.01*g/mole;
	G4Isotope* ih1 = new G4Isotope("Hydrogen",iz=1,n=1,a);
	a= 2.01*g/mole;
	G4Isotope* ih2 = new G4Isotope("Deuterium",iz=1,n=2,a);

	G4Element* elH = new G4Element(name="Hydrogen",symbol="H",2);
	elH->AddIsotope(ih1,.999);
	elH->AddIsotope(ih2,.001);

	a = 16.00*g/mole;
	G4Element* elO = new G4Element(name="Oxygen",symbol="O", z=8., a);

	density = 1.00*g/cm3;
	G4Material* Water = new G4Material(name="Water", density, nel=2);
	Water->AddElement(elH, 2);
	Water->AddElement(elO, 1);

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
// Add some properties to MPT  
//////////////////////////////////////////////////////////////////////////////

	myMPT.AddProperty("EFFIC", PPCKOV, EFFIC, NUMENTRIES);
	myMPT.AddProperty("ABSCO", PPCKOV, ABSCO, NUMENTRIES);
	myMPT.AddProperty("RINDEX", PPCKOV, RINDEX, NUMENTRIES);

//////////////////////////////////////////////////////////////////////////////
// test addition of MPT to G4Material 
//////////////////////////////////////////////////////////////////////////////

        G4cout << "\n Testing Addition of MPT to Water\n\n";

        // Add myMPT to Water  

        Water->SetMaterialPropertiesTable(&myMPT);

        // Get myMPT from the Water  

        G4MaterialPropertiesTable *anotherMPT;
        anotherMPT = Water->GetMaterialPropertiesTable();
	anotherMPT->DumpTable();
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
        while ( cin.get(ch) )
        {
                if (ch == '\n') break;
        }       
	G4cout << endl;
}

