// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsTableTest.cc,v 1.5 2001-03-05 12:00:27 gcosmo Exp $
// 
// ----------------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"

#include "g4std/iomanip"
 
int main()
{
//
// test the UnitsTable classes
//
// Build the Table of units
//
   G4UnitDefinition::BuildUnitsTable(); 
   G4UnitDefinition::PrintUnitsTable(); 

// Get internal value of an unit
//
   G4cout << "\n \t G4UnitDefinition::GetValueOf('Unit') \n";
      
   G4cout << " meter = " << G4UnitDefinition::GetValueOf("meter") << G4endl;
   G4cout << " cm    = " << G4UnitDefinition::GetValueOf("cm")    << G4endl;
   G4cout << " joule = " << G4UnitDefinition::GetValueOf("J")     << G4endl; 
   
// Get category of an unit
//
   G4cout << "\n \t G4UnitDefinition::GetCategory('Unit') \n";
      
   G4cout << " meter is " << G4UnitDefinition::GetCategory("m")    << G4endl;
   G4cout << " g     is " << G4UnitDefinition::GetCategory("gram") << G4endl;
   G4cout << " joule is " << G4UnitDefinition::GetCategory("J")    << G4endl; 
   G4cout << " ns    is " << G4UnitDefinition::GetCategory("ns")   << G4endl;
   
// Automatic conversion on output of a physical quantity
//
   G4cout << "\n \t G4BestUnit \n";
   G4cout.precision(3);   

   G4cout << " a = " << G4std::setw(4) << G4BestUnit (0.5*GeV ,"Energy") << G4endl;    
   G4cout << " b = " << G4std::setw(4) << G4BestUnit (0.15*MeV,"Energy") << G4endl;
   G4cout << " c = " << G4std::setw(4) << G4BestUnit (4000*MeV,"Energy") << G4endl;

   G4double x = -1000.*cm;   
   G4BestUnit d(x,"Length");         G4cout << " x = " << d << G4endl;
   G4BestUnit e(2e40*m,"Length");    G4cout << " e = " << e << G4endl;
   G4BestUnit f(DBL_MAX  ,"Energy"); G4cout << " f = " << f << G4endl;
   G4BestUnit h(0.,"Magnetic flux density"); G4cout << " h = " << h << G4endl;
   
   G4ThreeVector point(2*mm, 3*cm, 1*m);
   G4ThreeVector momen(3*MeV, 2*keV, 0.);
   G4cout << G4std::setw(6) << G4BestUnit (point, "Length") << G4endl;
   G4cout << G4std::setw(6) << G4BestUnit (momen, "Energy") << G4endl;
   
// Define new units
//
   new G4UnitDefinition("kg/m3","kg/m3","Volumic Mass",kg/m3);
   new G4UnitDefinition("g/cm3","g/cm3","Volumic Mass",g/cm3);

   G4double rho = 14.*mg/mm3;   
   G4cout << " rho = " << G4BestUnit (rho,"Volumic Mass") << G4endl;
   
   new G4UnitDefinition("g/cm2","g/cm2","Depth",g/cm2);
   
//   
//   G4UnitDefinition::PrintUnitsTable();           

   return 0;
}
