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
<<<<<<< HEAD
// $Id: RE05TrackerParametrisation.cc 66526 2012-12-19 13:41:33Z ihrivnac $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file RE05/src/RE05TrackerParametrisation.cc
/// \brief Implementation of the RE05TrackerParametrisation class
//

#include "RE05TrackerParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

RE05TrackerParametrisation::RE05TrackerParametrisation()
{

#include "RE05DetectorParameterDef.icc"

}

RE05TrackerParametrisation::~RE05TrackerParametrisation()
{;}

void RE05TrackerParametrisation::ComputeTransformation
(const G4int, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void RE05TrackerParametrisation::ComputeDimensions
(G4Tubs& trackerLayer, const G4int copyNo, const G4VPhysicalVolume*) const
{
  trackerLayer.SetInnerRadius(tracker_radius[copyNo]);
  trackerLayer.SetOuterRadius(tracker_radius[copyNo]+tracker_thick);
  trackerLayer.SetZHalfLength(tracker_length[copyNo]);
  trackerLayer.SetStartPhiAngle(trkTubs_sphi,false);
  trackerLayer.SetDeltaPhiAngle(trkTubs_dphi);
}
