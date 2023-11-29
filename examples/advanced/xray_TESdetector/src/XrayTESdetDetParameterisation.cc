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
/// \file DetParameterisation.cc
/// \brief Implementation of the DetParameterisation class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetDetParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
XrayTESdetDetParameterisation::XrayTESdetDetParameterisation(
  G4int NoChambers,
  G4double startX,           //  X of center of first
  G4double spacingX,         //  X spacing of centers
  G4double startY,           //  Y of center of first
  G4double spacingY          //  Y spacing of centers
){
  fNoChambers = NoChambers;
  fStartX = startX;
  fStartY = startY;
  fSpacingX = spacingX;
  fSpacingY = spacingY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetDetParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4double Xposition;
  G4double Yposition;
  G4int no_daug = physVol->GetLogicalVolume()->GetNoDaughters();

  for (G4int k=0; k<no_daug; k++)
  {
    physVol->GetLogicalVolume()->GetDaughter(k)->SetCopyNo(-copyNo);
  }
  std::ifstream input("pixelpos.txt");
  if ( !input )
  {
    G4cout << "Cannot open pixels positions' file\n"; exit( -1 );
  }

  std::vector<G4double> xpos; G4double b;
  std::vector<G4double> ypos; G4double c;

  while(!input.eof())
  {
    input >> b >> c;
    xpos.push_back(b);
    ypos.push_back(c);
  }

  input.close();

  Xposition = xpos[copyNo]*mm;
  Yposition = ypos[copyNo]*mm;

  G4ThreeVector origin(Xposition, Yposition, 0);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
