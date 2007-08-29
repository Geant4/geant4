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

#include "G4ScoringBoxParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

G4ScoringBoxParameterisation::G4ScoringBoxParameterisation(EAxis segaxis,
							   G4double motherDimension[3],
							   std::vector<G4double> & segpos)
  : fSegmentAxis(segaxis), fSegmentPositions(segpos) {

  for(int i = 0; i < 3; i++) fMotherDimensions[i] = motherDimension[i];
}

G4ScoringBoxParameterisation::~G4ScoringBoxParameterisation() {
  ;
}

void G4ScoringBoxParameterisation::ComputeTransformation(const G4int copyNo,
							G4VPhysicalVolume *physVol) const {
  G4ThreeVector trans;

  G4double pos;
  if(copyNo == 0)
    pos = fSegmentPositions[copyNo]/2.;
  else
    pos = (fSegmentPositions[copyNo-1] + fSegmentPositions[copyNo])/2.;

  if(fSegmentAxis == kXAxis) {
    trans[0] = pos - fMotherDimensions[0];
  }
  if(fSegmentAxis == kYAxis) {
    trans[1] = pos - fMotherDimensions[1];
  }
  if(fSegmentAxis == kZAxis) {
    trans[2] = pos - fMotherDimensions[2];
  }

  physVol->SetTranslation(trans);

}

void G4ScoringBoxParameterisation::ComputeDimensions(G4Box & meshElement,
						    const G4int copyNo,
						    const G4VPhysicalVolume*) const {

  G4double dims[3] = {fMotherDimensions[0],fMotherDimensions[1],fMotherDimensions[2]};

  G4double dim;
  if(copyNo == 0) 
    dim = fSegmentPositions[copyNo]/2.;
  else
    dim = (fSegmentPositions[copyNo] - fSegmentPositions[copyNo-1])/2.;

  if(fSegmentAxis == kXAxis) {
    dims[0] = dim;
  }
  if(fSegmentAxis == kYAxis) {
    dims[1] = dim;
  }
  if(fSegmentAxis == kZAxis) {
    dims[2] = dim;
  }

  meshElement.SetXHalfLength(dims[0]);
  meshElement.SetYHalfLength(dims[1]);
  meshElement.SetZHalfLength(dims[2]);
}
