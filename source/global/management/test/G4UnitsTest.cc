// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsTest.cc,v 1.4 1999-11-23 15:00:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

int main()
{
   G4cout << "          *** value of units *** " << G4endl << G4endl ;

   G4cout << "mm     = " << mm     << G4endl;

   G4cout << "cm     = " << cm     << G4endl;

   G4cout << "fermi  = " << fermi  << G4endl;
   
   G4cout << "cm2    = " << cm2    << G4endl;

   G4cout << "barn   = " << barn   << G4endl;

   G4cout << "second = " << s      << G4endl;

   G4cout << "joule  = " << joule  << G4endl;

   G4cout << "kg     = " << kg     << G4endl;

   G4cout << "watt   = " << watt   << G4endl;

   G4cout << "newton = " << newton << G4endl;

   G4cout << "pascal = " << pascal << G4endl;

   G4cout << "bar    = " << bar    << G4endl;

   G4cout << "coulomb= " << coulomb << G4endl;

   G4cout << "ampere = " << ampere << G4endl;

   G4cout << "volt   = " << volt   << G4endl;

   G4cout << "ohm    = " << ohm    << G4endl;

   G4cout << "farad  = " << farad  << G4endl;

   G4cout << "weber  = " << weber  << G4endl;

   G4cout << "tesla  = " << tesla  << G4endl;

   G4cout << "gauss  = " << gauss  << G4endl;

   G4cout << "henry  = " << henry  << G4endl;

   G4cout << "kelvin = " << kelvin << G4endl;

   G4cout << G4endl ;
//
// Physical Constants
//
   G4cout << "          *** Physical Constants *** " << G4endl << G4endl ;

   G4cout << "Avogadro    = " << Avogadro         << G4endl;

   G4cout << "c_light     = " << c_light          << G4endl;

   G4cout << "h_Planck    = " << h_Planck         << G4endl;

   G4cout << "hbar_Planck = " << hbar_Planck      << G4endl;

   G4cout << "hbarc       = " << hbarc            << G4endl;

   G4cout << "amu         = " << amu              << G4endl;
 
   G4cout << "mu0                     = " << mu0                      << G4endl;

   G4cout << "epsilon0                = " << epsilon0                 << G4endl;

   G4cout << "elm_coupling            = " << elm_coupling             << G4endl;

   G4cout << "fine_structure_const    = " << fine_structure_const     << G4endl;

   G4cout << "classic_electr_radius   = " << classic_electr_radius    << G4endl;

   G4cout << "electron_Compton_length = " << electron_Compton_length  << G4endl;

   G4cout << "Bohr_radius             = " << Bohr_radius              << G4endl;

   G4cout << "alpha_rcl2              = " << alpha_rcl2               << G4endl;

   G4cout << "twopi_mc2_rcl2          = " << twopi_mc2_rcl2           << G4endl;

   G4cout << "k_Boltzmann             = " << k_Boltzmann              << G4endl;

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
   
   G4cout << " meter = " << G4UnitDefinition::GetValueOf("meter") << G4endl;
   G4cout << " cm    = " << G4UnitDefinition::GetValueOf("cm")    << G4endl;
   G4cout << " joule = " << G4UnitDefinition::GetValueOf("J")     << G4endl; 
   
   G4double a = 0.5*GeV, b = 0.15*MeV, c = 4000*MeV;

   G4cout << "  a = " << G4BestUnit (a,"Energy") 
        << "  b = " << G4BestUnit (b,"Energy")
        << "  c = " << G4BestUnit (c,"Energy") << G4endl;

   return 0;
}
