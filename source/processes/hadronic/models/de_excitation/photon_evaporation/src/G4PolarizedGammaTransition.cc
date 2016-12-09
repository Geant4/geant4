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
// $Id: G4PolarizedGammaTransition.cc 85659 2014-11-03 10:59:10Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PolarizedGammaTransition
//
//      Author V.Ivanchenko 6 November 2015
//
// -------------------------------------------------------------------

#include "G4PolarizedGammaTransition.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationTransition.hh"

G4PolarizedGammaTransition::G4PolarizedGammaTransition() 
{
  fPolarization = new G4PolarizationTransition();
}

G4PolarizedGammaTransition::~G4PolarizedGammaTransition() 
{
  delete fPolarization;
}

void G4PolarizedGammaTransition::SampleDirection(
     G4Fragment*, G4double, G4int, G4int, G4int)
{
  fDirection = G4RandomDirection();
}
