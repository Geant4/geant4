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
//
// $Id: G4VPGeoDriver.hh,v 1.4 2002-08-29 15:30:50 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VPGeoDriver
//
// Class description:
//
// Used internally by importance sampling and scoring in a "parallel"
// geometry.
// It defines an interface to "drive" or "move" a particle in a 
// "parallel" geometry. 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VPGeoDriver_hh
#define G4VPGeoDriver_hh G4VPGeoDriver_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryCell.hh"

class G4VPGeoDriver
{

public:  // with description

  virtual ~G4VPGeoDriver(){}
  
  virtual G4GeometryCell LocateOnBoundary(const G4ThreeVector &aPosition, 
		                           const G4ThreeVector &aDirection) = 0;
    // The location of a track according to it's position
    // and direction in case the track crosses a boundary
    // of a "parallel" geometry.
    // Must be called in the PostDOIT of the ParallelTransportation.
    // (The track crosses the boundary if PostDOIT gets called.)
  
  virtual G4GeometryCell GetCurrentGeometryCell() const = 0;
    // get the current G4GeometryCell of the "parallel" geometry

  virtual G4double ComputeStepLengthInit(const G4ThreeVector &aPosition, 
                                         const G4ThreeVector &aDirection) = 0;
    // compute step length for a starting track. 

  virtual G4double ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
				                  const G4ThreeVector &aDirection) = 0;
    // compute the step length after a track crossed a boundary
    // in a "parallel" geometry
 
  virtual G4double ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
                                             const G4ThreeVector &aDirection) = 0;
    // compute step length when track moves inside a volume.  

};

#endif
