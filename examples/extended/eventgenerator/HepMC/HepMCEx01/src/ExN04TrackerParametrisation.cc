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
/// \file eventgenerator/HepMC/HepMCEx01/src/ExN04TrackerParametrisation.cc
/// \brief Implementation of the ExN04TrackerParametrisation class
//
//

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "ExN04TrackerParametrisation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04TrackerParametrisation::ExN04TrackerParametrisation()
 : G4VPVParameterisation()
{
#include "ExN04DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04TrackerParametrisation::~ExN04TrackerParametrisation()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04TrackerParametrisation::ComputeTransformation
                        (const G4int, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector origin;
  physVol-> SetTranslation(origin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN04TrackerParametrisation::ComputeDimensions
     (G4Tubs& trackerLayer, const G4int copyNo, const G4VPhysicalVolume*) const
{
  trackerLayer.SetInnerRadius(ftracker_radius[copyNo]);
  trackerLayer.SetOuterRadius(ftracker_radius[copyNo] + ftracker_thick);
  trackerLayer.SetZHalfLength(ftracker_length[copyNo]);
  trackerLayer.SetStartPhiAngle(ftrkTubs_sphi);
  trackerLayer.SetDeltaPhiAngle(ftrkTubs_dphi);
}
