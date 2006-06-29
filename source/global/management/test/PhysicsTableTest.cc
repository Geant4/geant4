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
// $Id: PhysicsTableTest.cc,v 1.6 2006-06-29 19:05:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//
// This program shows how to use the G4PhysicsTable and G4PhysicsVector. 
//

#include "G4ios.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include <cmath>

static const size_t PhysicsTableSize = 10;
static const size_t PhysicsVectorSize = 10;

int main()
{
 
  // Create an empty G4PhysicsTable Object
  G4PhysicsTable* pPhysTable;
  pPhysTable = new G4PhysicsTable();
 
  // reserve memory  
  pPhysTable->reserve(PhysicsTableSize);
  
  // create PhysicsVectors and store them into the PhysicsTable
  G4PhysicsVector* pVector;
  G4double e_min = 0.0;
  G4double e_max = 1.0;
  size_t   n_bin = PhysicsVectorSize;
  for (size_t i=0; i<PhysicsTableSize; i++){
    pVector = new G4PhysicsLinearVector(e_min, e_max, n_bin);
    pVector->PutComment("PhysicsLinearVector");

    // put values in PhysicsVector 
    //   value = std::sin(factor*energy)
    G4double factor = 0.1*G4double(i+1);
    for (size_t k=0; k<n_bin; k++){
      G4double eVal = pVector->GetLowEdgeEnergy(k);
      pVector->PutValue(k, std::sin(factor*eVal));
    }
    
    // put PhysicsVector in PhysicsTable
    pPhysTable->push_back(pVector);
  }
  
  // printout PhysicsTable
  // G4cout << *pPhysTable;

  // store PhysicsTable
  pPhysTable->StorePhysicsTable("testPhysTbleAscii.dat",true);

  // delete PhysicsTable
  pPhysTable->clearAndDestroy();
  delete pPhysTable;

  // create empty PhysicsTable 
  pPhysTable = new G4PhysicsTable();
 
  // retrieve PhysicsTable
  pPhysTable->RetrievePhysicsTable("testPhysTbleAscii.dat",true);

  // printout PhysicsTable 
  G4cout << *pPhysTable;
 
  // clear PhysicsTable
  pPhysTable->clearAndDestroy();
  delete pPhysTable;


  ///////////////////////////////////////////////////////////
  
  // create empty PhysicsTable 
  pPhysTable = new G4PhysicsTable();
  
  // create PhysicsVectors and store them into the PhysicsTable
  e_min =    0.1;
  e_max = 1.0e9;
  n_bin = PhysicsVectorSize;
  for (size_t j=0; j<PhysicsTableSize; j++){
    pVector = new G4PhysicsLogVector(e_min, e_max, n_bin);
    pVector->PutComment("PhysicsLogVector");

    // put values in PhysicsVector 
    //   value = std::log10(factor*energy)
    G4double factor = 0.1*G4double(j+1);
    for (size_t k=0; k<n_bin; k++){
      G4double eVal = pVector->GetLowEdgeEnergy(k);
      pVector->PutValue(k, std::log10(factor*eVal));
    }
    
    // put PhysicsVector in PhysicsTable
    pPhysTable->push_back(pVector);
  }
  
  // printout PhysicsTable
  // G4cout << *pPhysTable;

  // store PhysicsTable
  pPhysTable->StorePhysicsTable("testPhysTbleBin.dat");

  // delete PhysicsTable
  pPhysTable->clearAndDestroy();

  // retrieve PhysicsTable
   pPhysTable->RetrievePhysicsTable("testPhysTbleBin.dat");

  // printout PhysicsTable
  G4cout << *pPhysTable;
  

  return 0; 
} 

