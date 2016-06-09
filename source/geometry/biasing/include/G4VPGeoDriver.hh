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
//
// $Id: G4VPGeoDriver.hh,v 1.7 2006/06/29 18:16:45 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

class G4VPhysicalVolume;

class G4VPGeoDriver
{

public:  // with description

  G4VPGeoDriver();
  virtual ~G4VPGeoDriver();
  
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
