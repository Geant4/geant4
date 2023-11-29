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
//
// G4VSensitiveDetector
#include "G4SDChargedFilter.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector.
//
// Created: 2005-11-14  Tsukasa ASO.
//
///////////////////////////////////////////////////////////////////////////////

G4SDChargedFilter::G4SDChargedFilter(G4String name)
  : G4VSDFilter(name)
{}

G4bool G4SDChargedFilter::Accept(const G4Step* aStep) const
{
  return aStep->GetPreStepPoint()->GetCharge() != 0.;
}
