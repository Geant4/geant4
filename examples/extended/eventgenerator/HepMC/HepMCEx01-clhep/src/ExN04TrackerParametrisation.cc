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

#include "ExN04TrackerParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

ExN04TrackerParametrisation::ExN04TrackerParametrisation()
{

#include "ExN04DetectorParameterDef.icc"

}

ExN04TrackerParametrisation::~ExN04TrackerParametrisation()
{;}

void ExN04TrackerParametrisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void ExN04TrackerParametrisation::ComputeDimensions
(G4Tubs & trackerLayer, const G4int copyNo,
 const G4VPhysicalVolume * physVol) const
{
  trackerLayer.SetInnerRadius(tracker_radius[copyNo]);
  trackerLayer.SetOuterRadius(tracker_radius[copyNo]+tracker_thick);
  trackerLayer.SetZHalfLength(tracker_length[copyNo]);
  trackerLayer.SetStartPhiAngle(trkTubs_sphi);
  trackerLayer.SetDeltaPhiAngle(trkTubs_dphi);
}
