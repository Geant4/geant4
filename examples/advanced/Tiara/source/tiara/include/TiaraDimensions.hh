// $Id: TiaraDimensions.hh,v 1.2 2003-06-16 17:06:45 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
