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
//---------------------------------------------------------------------------
//
// ClassName:      G4HadronicParameters
//
// Author:         2018 Alberto Ribon
//
// Description:    Singleton to keep global hadronic parameters.
//
//                 For the time being, at least, offers only "getters" but
//                 not "setters", i.e. a recompilation is needed to change
//                 the default parameters.
//
// Modified:
//
//----------------------------------------------------------------------------

#include "G4HadronicParameters.hh"
#include <CLHEP/Units/PhysicalConstants.h>


G4HadronicParameters* G4HadronicParameters::sInstance = nullptr;


G4HadronicParameters* G4HadronicParameters::Instance() {
  if ( sInstance == nullptr ) {
    static G4HadronicParameters theHadronicParametersObject;
    sInstance = &theHadronicParametersObject;
  }
  return sInstance;
}


G4HadronicParameters::~G4HadronicParameters() {}


G4HadronicParameters::G4HadronicParameters() {
  fMaxEnergy = 100.0*CLHEP::TeV;
  fMinEnergyTransitionFTF_Cascade = 3.0*CLHEP::GeV;
  fMaxEnergyTransitionFTF_Cascade = 6.0*CLHEP::GeV;
  fMinEnergyTransitionQGS_FTF = 12.0*CLHEP::GeV;
  fMaxEnergyTransitionQGS_FTF = 25.0*CLHEP::GeV;
}

