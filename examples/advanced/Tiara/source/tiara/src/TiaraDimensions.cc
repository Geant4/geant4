// $Id: TiaraDimensions.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraDimensions.hh"

TiaraDimensions::TiaraDimensions() :
 
  // World volume
  worldHalfLength(330 * cm),
  worldHalfWidth(100 * cm),

  // basic dimensions
  targetPosZ(-325 * cm),
  distTargetWall(176 * cm),
  distTargetEndA(396 * cm),
  distTargetExperiment(401 * cm),
  distTargetEndB(516 * cm),
  
  // beam pipe
  pipeRadius((10.9/2) * cm),

  // radius iron A
  radiusIronA(26 * cm),
  
  // experiment width
  widthExperiment(120 * cm),

  // detector
  detectorRadius(12.7 / 2 * cm),
  detectorHalfHight(12.7 / 2 * cm),

  // source detector
  srcDetectorWidth(0.1 * cm)

{}

TiaraDimensions::~TiaraDimensions()
{}


