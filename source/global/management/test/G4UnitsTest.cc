// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsTest.cc,v 1.1 1999-01-07 16:09:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

int main()
{
   G4cout << "          *** value of units *** " << endl << endl ;

   G4cout << "mm     = " << mm     << endl;

   G4cout << "cm     = " << cm     << endl;

   G4cout << "fermi  = " << fermi  << endl;
   
   G4cout << "cm2    = " << cm2    << endl;

   G4cout << "barn   = " << barn   << endl;

   G4cout << "second = " << s      << endl;

   G4cout << "joule  = " << joule  << endl;

   G4cout << "kg     = " << kg     << endl;

   G4cout << "Watt   = " << Watt   << endl;

   G4cout << "newton = " << newton << endl;

   G4cout << "pascal = " << pascal << endl;

   G4cout << "bar    = " << bar    << endl;

   G4cout << "coulomb= " << coulomb << endl;

   G4cout << "ampere = " << ampere << endl;

   G4cout << "volt   = " << volt   << endl;

   G4cout << "ohm    = " << ohm    << endl;

   G4cout << "farad  = " << farad  << endl;

   G4cout << "weber  = " << weber  << endl;

   G4cout << "tesla  = " << tesla  << endl;

   G4cout << "gauss  = " << gauss  << endl;

   G4cout << "henry  = " << henry  << endl;

   G4cout << "kelvin = " << kelvin << endl;

   G4cout << endl ;
//
// Physical Constants
//
   G4cout << "          *** Physical Constants *** " << endl << endl ;

   G4cout << "Avogadro    = " << Avogadro         << endl;

   G4cout << "c_light     = " << c_light          << endl;

   G4cout << "h_Planck    = " << h_Planck         << endl;

   G4cout << "hbar_Planck = " << hbar_Planck      << endl;

   G4cout << "hbarc       = " << hbarc            << endl;

   G4cout << "amu         = " << amu              << endl;
 
   G4cout << "mu0                     = " << mu0                      << endl;

   G4cout << "epsilon0                = " << epsilon0                 << endl;

   G4cout << "elm_coupling            = " << elm_coupling             << endl;

   G4cout << "fine_structure_const    = " << fine_structure_const     << endl;

   G4cout << "classic_electr_radius   = " << classic_electr_radius    << endl;

   G4cout << "electron_Compton_length = " << electron_Compton_length  << endl;

   G4cout << "Bohr_radius             = " << Bohr_radius              << endl;

   G4cout << "alpha_rcl2              = " << alpha_rcl2               << endl;

   G4cout << "twopi_mc2_rcl2          = " << twopi_mc2_rcl2           << endl;

   G4cout << "k_Boltzmann             = " << k_Boltzmann              << endl;

//
// test the UnitsTable class
//
   new G4UnitDefinition(     "meter","m" ,"Length",m);
   new G4UnitDefinition("centimeter","cm","Length",cm);
   new G4UnitDefinition("millimeter","mm","Length",mm);
   new G4UnitDefinition("millimeter3","mm3","Volume",mm3);
   new G4UnitDefinition(    "electronvolt","eV" ,"Energy",eV);
   new G4UnitDefinition("kiloelectronvolt","keV","Energy",keV);
   new G4UnitDefinition("megaelectronvolt","MeV","Energy",MeV);
   new G4UnitDefinition("gigaelectronvolt","GeV","Energy",GeV);
   new G4UnitDefinition("joule"           ,"J"  ,"Energy",joule);
    
   G4UnitDefinition::PrintUnitsTable(); 
   
   G4cout << " meter = " << G4UnitDefinition::GetValueOf("meter") << endl;
   G4cout << " cm    = " << G4UnitDefinition::GetValueOf("cm")    << endl;
   G4cout << " joule = " << G4UnitDefinition::GetValueOf("J")     << endl; 
   
   G4double a = 0.5*GeV, b = 0.15*MeV, c = 4000*MeV;

   G4cout << "  a = " << G4BestUnit (a,"Energy") 
        << "  b = " << G4BestUnit (b,"Energy")
        << "  c = " << G4BestUnit (c,"Energy") << endl;

   return 0;
}
