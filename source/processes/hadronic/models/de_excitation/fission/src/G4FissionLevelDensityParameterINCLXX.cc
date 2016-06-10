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
// $Id: G4FissionLevelDensityParameterINCLXX.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by D. Mancusi (6th October 2014)
//

#include "G4FissionLevelDensityParameterINCLXX.hh"
#include "G4HadronicException.hh"


G4FissionLevelDensityParameterINCLXX::G4FissionLevelDensityParameterINCLXX() :
  afanLow(1.02),
  afanHigh(1.02),
  ZLow(84),
  ZHigh(89)
{
  UpdateAfanSlope();
}

G4FissionLevelDensityParameterINCLXX::~G4FissionLevelDensityParameterINCLXX()
{}


G4double G4FissionLevelDensityParameterINCLXX::
LevelDensityParameter(G4int A, G4int Z, G4double U) const 
{
  G4double EvapLDP = 
    theEvaporationLevelDensityParameter.LevelDensityParameter(A,Z,U);

  if(Z >= ZHigh)     { EvapLDP *= afanHigh; }
  else if(Z <= ZLow) { EvapLDP *= afanLow; }
  else               { EvapLDP *= (afanLow + afanSlope*(Z-ZLow)); }

  return EvapLDP;

}

void G4FissionLevelDensityParameterINCLXX::UpdateAfanSlope() {
  if(ZHigh!=ZLow)
    afanSlope = (afanHigh-afanLow)/((double)(ZHigh-ZLow));
  else
    afanSlope = 0.;
}
