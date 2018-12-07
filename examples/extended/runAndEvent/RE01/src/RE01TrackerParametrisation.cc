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
/// \file runAndEvent/RE01/src/RE01TrackerParametrisation.cc
/// \brief Implementation of the RE01TrackerParametrisation class
//
//

#include "RE01TrackerParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01TrackerParametrisation::RE01TrackerParametrisation()
  : G4VPVParameterisation()
{
#include "RE01DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01TrackerParametrisation::~RE01TrackerParametrisation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01TrackerParametrisation
::ComputeTransformation(const G4int, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void RE01TrackerParametrisation
::ComputeDimensions(G4Tubs& trackerLayer, const G4int copyNo, 
                    const G4VPhysicalVolume*) const
{
  trackerLayer.SetInnerRadius(fTracker_radius[copyNo]);
  trackerLayer.SetOuterRadius(fTracker_radius[copyNo]+fTracker_thick);
  trackerLayer.SetZHalfLength(fTracker_length[copyNo]);
  trackerLayer.SetStartPhiAngle(fTrkTubs_sphi);
  trackerLayer.SetDeltaPhiAngle(fTrkTubs_dphi);
}
