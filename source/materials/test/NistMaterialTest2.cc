//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: NistMaterialTest2.cc,v 1.2 2005-02-22 10:11:09 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//
//
//  Test the construction of materials from the NIST data base
//
#include "G4NistManager.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

int main() {

G4int z;
G4bool buildIsotopes;

G4UnitDefinition::BuildUnitsTable();

// initialise NIST data base
//
G4NistManager*  nistMat = G4NistManager::Instance();

// define Elements
//
nistMat->FindOrBuildElement (z=1);
nistMat->FindOrBuildElement (z=6);
nistMat->FindOrBuildElement (z=7, buildIsotopes=false);
nistMat->FindOrBuildElement ("O");
nistMat->FindOrBuildElement ("Si", buildIsotopes=false);
nistMat->FindOrBuildElement ("Fe");
nistMat->FindOrBuildElement ("U");

G4cout << *(G4Isotope::GetIsotopeTable()) << G4endl;
G4cout << *(G4Element::GetElementTable()) << G4endl;

 
// define Materials
//
G4Material* Al =
nistMat->FindOrBuildMaterial ("G4_Al");
nistMat->FindOrBuildMaterial ("G4_lAr", buildIsotopes=false);
nistMat->FindOrBuildMaterial ("G4_Cu" , buildIsotopes=false);
G4Material* Pb =
nistMat->FindOrBuildMaterial ("G4_Pb");

G4Material* H2O =
nistMat->FindOrBuildMaterial ("G4_WATER");
nistMat->FindOrBuildMaterial ("G4_POLYSTYRENE");
nistMat->FindOrBuildMaterial ("G4_SILICON_DIOXIDE");
nistMat->FindOrBuildMaterial ("G4_AIR", buildIsotopes=false);
// build new materials from scratch
//
G4String gelSiElm[3]  = {"H", "O", "Si"};
G4int    gelSiAtom[3] = { 4,   4,   1};
G4double gelSiW[3]    = {0.0416, 0.6661, 0.2923};

std::vector<G4String> elm;
std::vector<G4int> atom;
std::vector<G4double> w;
for (G4int j=0; j<3; j++) { 
   elm.push_back(gelSiElm[j]); atom.push_back(gelSiAtom[j]);
   w.push_back(gelSiW[j]);
}

// by atom count
G4Material* Aerog =    
nistMat->ConstructNewMaterial ("gel_silicate1", elm, atom, 0.2);

// by mass fraction
nistMat->ConstructNewMaterial ("gel_silicate2", elm, w,    0.2);

// print G4MaterialTable
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

// print additional informations
//
G4cout << " Nuclear interaction length of Aluminium: " 
       << Al->GetNuclearInterLength()/cm << " cm" << G4endl;
       
G4cout << " Nuclear interaction length of Lead: " 
       << Pb->GetNuclearInterLength()/cm << " cm" << G4endl;
              
G4cout << " Nuclear interaction length of Water: " 
       << H2O->GetNuclearInterLength()/cm << " cm" << G4endl;
       
G4cout << " Nuclear interaction length of Aerogel: " 
       << Aerog->GetNuclearInterLength()/cm << " cm" << G4endl;
                           
return EXIT_SUCCESS;
}
