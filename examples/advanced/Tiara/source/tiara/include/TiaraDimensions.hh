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
// $Id: TiaraDimensions.hh,v 1.4 2006/06/29 15:43:34 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
//
// Class TiaraDimensions
//

#ifndef TiaraDimensions_hh
#define TiaraDimensions_hh TiaraGeometry_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

class TiaraDimensions {
public:
  TiaraDimensions();
  ~TiaraDimensions();

  // World volume
  G4double worldHalfLength;
  G4double worldHalfWidth;

  // basic dimensions
  G4double targetPosZ;
  G4double distTargetWall;
  G4double distTargetEndA;
  G4double distTargetExperiment;
  G4double distTargetEndB;
  
  // beam pipe
  G4double pipeRadius;

  // radius iron A
  G4double radiusIronA;
  
  // experiment width
  G4double widthExperiment;

  // detector
  G4double detectorRadius;
  G4double detectorHalfHight;

  // source detector
  G4double srcDetectorWidth;

};

#endif
