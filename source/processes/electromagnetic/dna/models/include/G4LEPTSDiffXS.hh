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
#ifndef G4LEPTSDiffXS_h
#define G4LEPTSDiffXS_h 1

#include <string>
#include "G4Types.hh"

class G4LEPTSDiffXS {

public:

  G4LEPTSDiffXS(std::string);   // Constructor

  void readDXS();    // Read file
  void BuildCDXS();
  void BuildCDXS(G4double, G4double);
  void NormalizeCDXS();
  void InterpolateCDXS();
  void PrintDXS(int);

  G4double SampleAngle(G4double);
  G4double SampleAngleMT(G4double, G4double);
  G4double SampleAngleEthylene(G4double, G4double);
  G4bool IsFileFound() const {
    return bFileFound;
  }

private:
  std::string fileName;
  int NumAng;
  int INumAng;
  int NumEn;
  char DXSTypeName[8];
  int DXSType;
  G4double Eb[100];
  //  G4double DXS[100][190], CDXS[100][190], IDXS[100][19000], ICDXS[100][19000];
  G4double DXS[100][190], CDXS[100][190], ICDXS[100][19000];
  //  G4double KT[100][190],  CKT[100][190],  IKT[100][19000];
  G4double KT[100][190],  IKT[100][19000];

  G4bool bFileFound;
};

#endif
