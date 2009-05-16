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
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4
// test for EnergyLoss Range InverseRange Lambda and other Tables, 
// A. Bagulya 14 May 2009
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc,char** argv)
{
  // Control on input

  if(argc < 2) {
    G4cout << "Input parameters are not specified! Exit" << G4endl;
    exit(1);
  }


  std::string path = "/home/bagoulia/geant4/examples/extended/electromagnetic/TestEm3/";
  std::string dir1 = "physdata_100bins/";
  std::string ff   = argv[1];
  std::string fName1 = path + dir1 + ff;

   std::ifstream in1;
   in1.open(fName1.c_str());
   if( !in1.is_open()) {
     G4cout << "Input file<" << fName1 << "> does not exist! Exit" << G4endl;
     exit(1);
   }

   G4PhysicsTable* t1 = new G4PhysicsTable();
   t1->RetrievePhysicsTable(fName1, true);

   G4PhysicsVector* V1 = (*t1)(1);
   
   //   V1->SetSpline(true);
   std::string d_st = "physdata_";
   std::string num_bins = argv[2];
   G4int n_spl = atoi(argv[3]);
   std::string s_spl = "spl";
   std::string d_m;
   if (n_spl) d_m = "bins_spl/";
   else d_m = "bins/";

   std::string dir2 = d_st + num_bins + d_m;
   std::string fName2 = path + dir2 + ff;

   std::ifstream in2;
   in2.open(fName2.c_str());

   if( !in2.is_open()) {
     G4cout << "Input file<" << fName2 << "> does not exist! Exit" << G4endl;
     exit(1);
   }

   G4double n = V1->GetVectorLength();
   G4double diff_max = 0.0;
   G4double e, y1, y2, diff;

   G4cout << "n = " << n <<G4endl;

   G4PhysicsTable* t2 = new G4PhysicsTable();
   t2->RetrievePhysicsTable(fName2, true);
 
   G4PhysicsVector* V2 = (*t2)(1); 

   std::string d_out = "data/";
   std::string typeFile = ".out";
   std::string c = "_";
   std::string bi = "bins";
 
   std::string asciiFileName;
   if (n_spl) {
     V2->SetSpline(true);
     asciiFileName = d_out + ff + c + num_bins + bi + c + s_spl + typeFile;
   } else {
     asciiFileName = d_out + ff + c + num_bins + bi + typeFile;
   }

   std::ofstream asciiFile;
   asciiFile.open(asciiFileName.c_str(), std::ios::out);
   if(asciiFile.is_open()) {
     asciiFile << " Energy(Mev) ||    Y1    ||    Y2      ||  Diff " << G4endl;
   } else {
     G4cout << "ERROR file <" << asciiFileName << "> is not opened" << G4endl;
     exit(1);
   } 

     G4bool b;
   for (G4int i = 0; i < n; i++) {
     e = V1->GetLowEdgeEnergy(i); 
     y1 = (*V1)[i];
     y2 = V2->GetValue(e, b);
     diff = std::fabs((1.0 - y2/y1)*100);
     if (diff > diff_max) diff_max = diff;

      asciiFile << std::setprecision(5)
                << std::setiosflags(std::ios::right);
      asciiFile << e;
      asciiFile << "    ";
      asciiFile << std::setiosflags(std::ios::right);
      asciiFile << y1;
      asciiFile << "    ";
      asciiFile << std::setiosflags(std::ios::right);
      asciiFile << y2;
      asciiFile << "    ";
      asciiFile << std::setiosflags(std::ios::right);
      asciiFile << diff
                << G4endl;

   }
   G4cout << "diff_max = " << diff_max << G4endl;

    delete V1;
    delete V2;
    exit(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

