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

















