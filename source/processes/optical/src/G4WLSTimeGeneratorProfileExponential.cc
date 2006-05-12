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
//
// GEANT4 Class file
//
//
// File name:     G4WLSTimeGeneratorProfileExponential.cc
//
// Author:        Pedro Rodrigues, Andreia Trindade
// 
// Creation date: 2006-05-07
//
// Modifications:
//
// Class Description: 
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4WLSTimeGeneratorProfileExponential.hh"
#include "Randomize.hh"
//    

G4WLSTimeGeneratorProfileExponential::G4WLSTimeGeneratorProfileExponential(const G4String& name):G4VWLSTimeGeneratorProfile(name)
{;}

//    

G4WLSTimeGeneratorProfileExponential::~G4WLSTimeGeneratorProfileExponential() 
{;}

//

G4double G4WLSTimeGeneratorProfileExponential::GenerateTime(const G4double time_constant)
{
  G4double time = 0;
  time = -std::log(G4UniformRand())*time_constant;
  return time;
}

G4double G4WLSTimeGeneratorProfileExponential::GenerateTime(const G4MaterialPropertiesTable*){
  // This method is not currently in use
  return 0;
}
