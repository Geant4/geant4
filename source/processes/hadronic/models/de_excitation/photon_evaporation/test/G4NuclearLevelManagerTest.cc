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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4NuclearLevelManagerTest.cc 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 27 October 1998
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <assert.h>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevel.hh"
#include "G4DataVector.hh"

int main()
{

  G4int Z;
  G4int A;

  G4cout << "Enter Z and A" << G4endl;
  G4cin >> Z >> A;

  assert (Z > 0);
  assert (A > 0);
  assert (A > Z);

  G4int iter;
  G4cout << "Enter number of iterations " << G4endl;
  G4cin >> iter;

   G4int i;
   for (i=0; i<iter; i++)
     {

       // NuclearLevelManager for this (Z,A) material
       //       G4NuclearLevelManager levelManager(Z,A);
       //       G4NuclearLevelManager* levelManager = new G4NuclearLevelManager(Z,A);
       G4NuclearLevelManager* levelManager = new G4NuclearLevelManager();
       levelManager->SetNucleus(Z,A);
       //       G4cout << i <<") ---- G4NuclearLevelManager created ----" << G4endl;
       
       // Is it a valid isotope?
       if (! (levelManager->IsValid(Z,A)))
	   	 { G4cout << "This (Z,A) has no ENSDF levels and gammas " << G4endl; }

       levelManager->PrintAll();
       delete levelManager;
     }

  return EXIT_SUCCESS;
}

















