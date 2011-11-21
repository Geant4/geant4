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
// $Id: test11_hadronic_exerciser.cc,v 1.2 2009-10-20 07:36:13 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 080901 Add Thermal materials and hydogen from 1H isotope in test.
//

#include "globals.hh"
#include "G4UnitsTable.hh"
#include <vector>

static void OutputCases
(G4int N,  // Number of events per case.
 const std::vector <G4String> & particleNameList,
 const std::vector <G4double> & energyList,
 const std::vector <G4String> & materialNameList) {

  for (size_t iMaterial = 0;
       iMaterial < materialNameList.size ();
       iMaterial++) {
    for (size_t iEnergy = 0;
	 iEnergy < energyList.size ();
	 iEnergy++) {
      for (size_t iParticle = 0;
	   iParticle < particleNameList.size ();
	   iParticle++) {

	G4cout
	  << "\n#\n# " << particleNameList [iParticle]
	  << " at " << G4BestUnit (energyList [iEnergy], "Energy")
	  << " in " << materialNameList [iMaterial]
	  << "\n#";

	G4cout
	  << "\n/gun/particle " << particleNameList [iParticle]
	  << "\n/gun/energy " <<  G4BestUnit (energyList [iEnergy], "Energy")
	  << "\n/mydet/SelectMaterial " << materialNameList [iMaterial]
 	  << "\n/run/beamOn " << N;

      }
    }
  }
}


int main (int argc, char** argv) {

  G4int N = 1;
  if (argc > 1) {
    if (strcmp (argv[1], "large_N") == 0) {
      N = 20;
    }
  }

  G4UnitDefinition::BuildUnitsTable();

  G4cout <<
    "#"
    "\n# Auto-generated test input file for test11 hadronics."
    "\n#"
    "\n/control/verbose 2"
    "\n# /run/verbose 2"
    "\n/run/setCut 1 m"
    "\n/run/initialize"
    "\n/gun/direction 0 0 1";

  std::vector <G4String> particleNameList;
  particleNameList.push_back ("proton");
  particleNameList.push_back ("neutron");
  particleNameList.push_back ("pi+");
  particleNameList.push_back ("pi-");
  particleNameList.push_back ("kaon+");
  particleNameList.push_back ("kaon-");
  particleNameList.push_back ("kaon0S");
  particleNameList.push_back ("kaon0L");

  std::vector <G4double> energyList;
  energyList.push_back (10 * GeV);

  std::vector <G4String> materialNameList;
  materialNameList.push_back ("Pb");
  materialNameList.push_back ("Al");
  materialNameList.push_back ("Air");

  OutputCases (N, particleNameList, energyList, materialNameList);

// Below T. K. add  for neutron HP
  std::vector <G4String> TK_particleNameList;
  TK_particleNameList.push_back ("neutron");
  std::vector <G4double> TK_energyList;
  TK_energyList.push_back (   10 * MeV);
  TK_energyList.push_back (    1 * MeV);
  TK_energyList.push_back (  100 * keV);
  TK_energyList.push_back (   10 * keV);
  TK_energyList.push_back (    1 * keV);
  TK_energyList.push_back (  100 *  eV);
  TK_energyList.push_back (   10 *  eV);
  TK_energyList.push_back (    1 *  eV);
  TK_energyList.push_back (  0.1 *  eV);
  TK_energyList.push_back ( 0.01 *  eV);
  TK_energyList.push_back (0.001 *  eV);

  std::vector <G4String> TK_materialNameList;
   TK_materialNameList.push_back ( "G4_H" );
/*
   TK_materialNameList.push_back ( "G4_He" );
   TK_materialNameList.push_back ( "G4_Li" );
   TK_materialNameList.push_back ( "G4_Be" );
   TK_materialNameList.push_back ( "G4_B" );
*/
   TK_materialNameList.push_back ( "G4_C" );
/*
   TK_materialNameList.push_back ( "G4_N" );
   TK_materialNameList.push_back ( "G4_O" );
   TK_materialNameList.push_back ( "G4_F" );
   TK_materialNameList.push_back ( "G4_Ne" );
   TK_materialNameList.push_back ( "G4_Na" );
   TK_materialNameList.push_back ( "G4_Mg" );
*/
   TK_materialNameList.push_back ( "G4_Al" );
/*
   TK_materialNameList.push_back ( "G4_Si" );
   TK_materialNameList.push_back ( "G4_P" );
   TK_materialNameList.push_back ( "G4_S" );
   TK_materialNameList.push_back ( "G4_Cl" );
   TK_materialNameList.push_back ( "G4_Ar" );
   TK_materialNameList.push_back ( "G4_K" );
   TK_materialNameList.push_back ( "G4_Ca" );
   TK_materialNameList.push_back ( "G4_Sc" );
   TK_materialNameList.push_back ( "G4_Ti" );
   TK_materialNameList.push_back ( "G4_V" );
   TK_materialNameList.push_back ( "G4_Cr" );
   TK_materialNameList.push_back ( "G4_Mn" );
*/
   TK_materialNameList.push_back ( "G4_Fe" );
/*
   TK_materialNameList.push_back ( "G4_Co" );
   TK_materialNameList.push_back ( "G4_Ni" );
   TK_materialNameList.push_back ( "G4_Cu" );
   TK_materialNameList.push_back ( "G4_Zn" );
   TK_materialNameList.push_back ( "G4_Ga" );
   TK_materialNameList.push_back ( "G4_Ge" );
   TK_materialNameList.push_back ( "G4_As" );
   TK_materialNameList.push_back ( "G4_Se" );
   TK_materialNameList.push_back ( "G4_Br" );
   TK_materialNameList.push_back ( "G4_Kr" );
   TK_materialNameList.push_back ( "G4_Rb" );
   TK_materialNameList.push_back ( "G4_Sr" );
   TK_materialNameList.push_back ( "G4_Y" );
   TK_materialNameList.push_back ( "G4_Zr" );
   TK_materialNameList.push_back ( "G4_Nb" );
   TK_materialNameList.push_back ( "G4_Mo" );
   TK_materialNameList.push_back ( "G4_Tc" );
   TK_materialNameList.push_back ( "G4_Ru" );
   TK_materialNameList.push_back ( "G4_Rh" );
   TK_materialNameList.push_back ( "G4_Pd" );
*/
   TK_materialNameList.push_back ( "G4_Ag" );
/*
   TK_materialNameList.push_back ( "G4_Cd" );
   TK_materialNameList.push_back ( "G4_In" );
   TK_materialNameList.push_back ( "G4_Sn" );
   TK_materialNameList.push_back ( "G4_Sb" );
   TK_materialNameList.push_back ( "G4_Te" );
   TK_materialNameList.push_back ( "G4_I" );
   TK_materialNameList.push_back ( "G4_Xe" );
   TK_materialNameList.push_back ( "G4_Cs" );
   TK_materialNameList.push_back ( "G4_Ba" );
   TK_materialNameList.push_back ( "G4_La" );
   TK_materialNameList.push_back ( "G4_Ce" );
   TK_materialNameList.push_back ( "G4_Pr" );
   TK_materialNameList.push_back ( "G4_Nd" );
   TK_materialNameList.push_back ( "G4_Pm" );
   TK_materialNameList.push_back ( "G4_Sm" );
   TK_materialNameList.push_back ( "G4_Eu" );
   TK_materialNameList.push_back ( "G4_Gd" );
   TK_materialNameList.push_back ( "G4_Tb" );
   TK_materialNameList.push_back ( "G4_Dy" );
   TK_materialNameList.push_back ( "G4_Ho" );
   TK_materialNameList.push_back ( "G4_Er" );
   TK_materialNameList.push_back ( "G4_Tm" );
   TK_materialNameList.push_back ( "G4_Yb" );
   TK_materialNameList.push_back ( "G4_Lu" );
   TK_materialNameList.push_back ( "G4_Hf" );
   TK_materialNameList.push_back ( "G4_Ta" );
   TK_materialNameList.push_back ( "G4_W" );
   TK_materialNameList.push_back ( "G4_Re" );
   TK_materialNameList.push_back ( "G4_Os" );
   TK_materialNameList.push_back ( "G4_Ir" );
   TK_materialNameList.push_back ( "G4_Pt" );
   TK_materialNameList.push_back ( "G4_Au" );
   TK_materialNameList.push_back ( "G4_Hg" ); // No data HP
   TK_materialNameList.push_back ( "G4_Tl" );
*/
   TK_materialNameList.push_back ( "G4_Pb" );
/*
   TK_materialNameList.push_back ( "G4_Bi" );
   TK_materialNameList.push_back ( "G4_Po" );
   //TK_materialNameList.push_back ( "G4_At" );
   //TK_materialNameList.push_back ( "G4_Rn" ); // No data HP
   //TK_materialNameList.push_back ( "G4_Fr" ); // No data HP
   //TK_materialNameList.push_back ( "G4_Ra" );
   //TK_materialNameList.push_back ( "G4_Ac" );
   //TK_materialNameList.push_back ( "G4_Th" );
   //TK_materialNameList.push_back ( "G4_Pa" );
*/
   TK_materialNameList.push_back ( "G4_U" );
   //TK_materialNameList.push_back ( "G4_Np" );
   //TK_materialNameList.push_back ( "G4_Pu" );
   //TK_materialNameList.push_back ( "G4_Am" );
   //TK_materialNameList.push_back ( "G4_Cm" );
   //TK_materialNameList.push_back ( "G4_Bk" );
   //TK_materialNameList.push_back ( "G4_Cf" ); // No data HP even with Am

//080901
   TK_materialNameList.push_back ( "Hydrogen1" );

//111025 Add Nist Materials
   TK_materialNameList.push_back ( "Water_TS" );
   TK_materialNameList.push_back ( "G4_WATER" ); 
   TK_materialNameList.push_back ( "Polyethylene_TS" );
   TK_materialNameList.push_back ( "G4_POLYETHYLENE" ); 
   TK_materialNameList.push_back ( "Graphite_TS" );
   TK_materialNameList.push_back ( "G4_GRAPHITE" ); 

   TK_materialNameList.push_back ( "Heavy_Water_TS" ); 
   //TK_materialNameList.push_back ( "ZrH_TS" ); 
   TK_materialNameList.push_back ( "Beryllium_Oxide_TS" ); 
   TK_materialNameList.push_back ( "G4_BERYLLIUM_OXIDE" ); 
   //TK_materialNameList.push_back ( "Uranium_Dioxide_TS" ); 
   //TK_materialNameList.push_back ( "G4_URANIUM_OXIDE" ); 

   TK_materialNameList.push_back ( "Be_TS" ); 
   TK_materialNameList.push_back ( "Al_TS" ); 
   TK_materialNameList.push_back ( "Fe_TS" ); 

//111101 Add Excited Isomer Materials
   TK_materialNameList.push_back ( "Water_TS" );
   TK_materialNameList.push_back ( "27Co58m1" );
   TK_materialNameList.push_back ( "47Ag110m1" );
   TK_materialNameList.push_back ( "48Cd115m1" );
   TK_materialNameList.push_back ( "52Te127m1" );
   TK_materialNameList.push_back ( "52Te129m1" );
   TK_materialNameList.push_back ( "61Pm148m1" );
   TK_materialNameList.push_back ( "67Ho166m1" );
   //TK_materialNameList.push_back ( "95Am242m1" );
   //TK_materialNameList.push_back ( "95Am244m1" );

/*
  if (argc > 1) {
    if (strcmp (argv[1], "large_N") == 0) {
      N = 10;
    }
  }
*/
  OutputCases (N, TK_particleNameList, TK_energyList, TK_materialNameList);

  G4cout << G4endl;
}
