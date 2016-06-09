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
// $Id: TiaraDimensions.cc,v 1.3 2003/06/25 09:12:59 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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


