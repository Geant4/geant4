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
// Anton.Lechner@cern.ch
//
// Environment variable TESTTARGET must be specified (see GNUmakefile)
//
//
// The following functions of the class G4ExtDEDXTable are tested:
//       IsApplicable()
//       AddPhysicsVector()
//       RemovePhysicsVector()
//       GetPhysicsVector()
//       StorePhysicsTable()
//       RetrievePhysicsTable()
//
// Note: This test creates a file "usertable.dat". You may want to delete it
//       after running the test.

#include "globals.hh"
#include "G4ExtDEDXTable.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"


void Test(const G4String& fileName) {

  // ########################################################################
  // # Creating physics vectors
  // ########################################################################

  G4double factor = MeV * cm2 / (0.001 * g);

  G4double E[31] = {.025,.03,.04,.05,.06,.07,.08,.09,.1,.15,.2,.25,.3,
                    .4,.5,.6,.7,.8,.9,1,1.5,2,2.5,3,4,5,6,7,8,9,10};

  // Lithium ions in hydrogen
  G4double H_3[31]={8.31,8.919,9.838,10.46,10.87,11.14,11.3,11.4,11.43,
                    11.21,10.73,10.21,9.696,8.773,7.989,7.326,6.761,6.275,
                    5.853,5.485,4.175,3.377,2.843,2.46,1.948,1.62,1.393,
                    1.225,1.095,0.9921,0.9082};

  G4LPhysicsFreeVector* physVecZ3MatH = 
                 new G4LPhysicsFreeVector(31, E[0] * MeV, E[30] * MeV);

  for(size_t i = 0; i < 31; i++) {
    physVecZ3MatH -> PutValues(i, E[i] * MeV, H_3[i] * factor);
  }
  physVecZ3MatH -> SetSpline(true);

  // Lithium ions in tissue
  G4double A_3[31]={2.748,2.992,3.386,3.674,3.877,4.016,4.108,4.166,4.2,
                    4.186,4.059,3.906,3.751,3.458,3.2,2.974,2.777,2.603,
                    2.449,2.313,1.807,1.485,1.263,1.1,0.8801,0.7372,0.6368,
                    0.5623,0.5046,0.4586,0.4209}; 

  G4LPhysicsFreeVector* physVecZ3MatTissue = 
                 new G4LPhysicsFreeVector(31, E[0] * MeV, E[30] * MeV);

  for(size_t i = 0; i < 31; i++) {
      physVecZ3MatTissue -> PutValues(i, E[i] * MeV, A_3[i] * factor);
  }
  physVecZ3MatTissue -> SetSpline(true);

  // ########################################################################
  // # Test 1:
  // ########################################################################

  G4ExtDEDXTable* tableA = new G4ExtDEDXTable();

  G4cout << G4endl << "##### Testing AddPhysicsVector() #####" 
         << G4endl;

  G4cout << "***Stopping power of hydrogen for Li ions***" << G4endl;

  G4cout << "  Vector not available. Is applicable: "; 
  if( tableA -> IsApplicable(3, "G4_H") ) 
     G4cout << "Yes. Test failed" << G4endl; 
  else 
     G4cout << "No. Test passed" << G4endl; 

  G4cout << "  Adding vector. Success: "; 
  if( tableA -> AddPhysicsVector(physVecZ3MatH, 3, "G4_H", 1) )
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "  Vector added. Is applicable: "; 
  if( tableA -> IsApplicable(3, "G4_H") ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "  Vector added. Is applicable: "; 
  if( tableA -> IsApplicable(3, 1) ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "***Stopping power of tissue for Li ions***" << G4endl;
  
  G4cout << "  Vector not available. Is applicable: "; 
  if( tableA -> IsApplicable(3, "G4_A-150_TISSUE") ) 
     G4cout << "Yes. Test failed" << G4endl; 
  else 
     G4cout << "No. Test passed" << G4endl; 

  G4cout << "  Adding vector. Success: "; 
  if( tableA -> AddPhysicsVector(physVecZ3MatTissue, 3, "G4_A-150_TISSUE") )
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "  Vector added. Is applicable: "; 
  if( tableA -> IsApplicable(3, "G4_A-150_TISSUE") ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  // ########################################################################
  // # Test 2:
  // ########################################################################

  G4cout << G4endl << "##### Testing GetPhysicsVector() #####" 
         << G4endl;

  G4PhysicsVector* vector = 0;

  vector = tableA -> GetPhysicsVector(3, "G4_H");

  G4cout << "  Vector retrieved. Found: "; 
  if( physVecZ3MatH == vector ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  vector = tableA -> GetPhysicsVector(3, 1);

  G4cout << "  Vector retrieved. Found: ";
  if( physVecZ3MatH == vector) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  vector = tableA -> GetPhysicsVector(3, "G4_A-150_TISSUE");
 
  G4cout << "  Vector retrieved. Found: ";
  if( physVecZ3MatTissue == vector) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 


  // ########################################################################
  // # Test 3:
  // ########################################################################

  G4cout << G4endl << "##### Testing StorePhysicsTable() #####" 
         << G4endl;

  G4cout << "  Storing table: ";
  if( tableA -> StorePhysicsTable(fileName) ) 
     G4cout << " Test passed" << G4endl;
  else 
     G4cout << " Test failed" << G4endl; 

  // ########################################################################
  // # Test 4:
  // ########################################################################

  G4cout << G4endl << "##### Testing RemovePhysicsVector() #####" 
         << G4endl;

  G4cout << "  Removing vector. Success "; 
  if( tableA -> RemovePhysicsVector(3, "G4_H") )
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "  Vector was removed. Is applicable: "; 
  if( tableA -> IsApplicable(3, "G4_H") ) 
     G4cout << "Yes. Test failed" << G4endl; 
  else 
     G4cout << "No. Test passed" << G4endl; 

  G4cout << "  Vector was removed. Is applicable: "; 
  if( tableA -> IsApplicable(3, 1) ) 
     G4cout << "Yes. Test failed" << G4endl; 
  else 
     G4cout << "No. Test passed" << G4endl; 

  delete tableA;

  // ########################################################################
  // # Test 5:
  // ########################################################################

  G4ExtDEDXTable* tableB = new G4ExtDEDXTable();

  G4cout << G4endl << "##### Testing RetrievePhysicsTable() #####" 
         << G4endl;

  G4cout << "  Vector not available. Is applicable: "; 
  if( tableB -> IsApplicable(3, "G4_H") ) 
     G4cout << "Yes. Test failed" << G4endl; 
  else 
     G4cout << "No. Test passed" << G4endl; 

  G4cout << "  Retrieving table: ";
  if( tableB -> RetrievePhysicsTable(fileName) ) 
     G4cout << "Test passed" << G4endl;
  else 
     G4cout << "Test failed" << G4endl;

  G4cout << "  Vector added. Is applicable: "; 
  if( tableB -> IsApplicable(3, "G4_H") ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << "  Vector added. Is applicable: "; 
  if( tableB -> IsApplicable(3, 1) ) 
     G4cout << "Yes. Test passed" << G4endl; 
  else 
     G4cout << "No. Test failed" << G4endl; 

  G4cout << G4endl;

  delete tableB;
}


int main() {  


  Test("usertable.dat");

  return EXIT_SUCCESS;
} 

