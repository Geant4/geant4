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
// $Id: G4VPGeoDriver.hh,v 1.2 2002-04-09 16:23:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VPGeoDriver
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VPGeoDriver_hh
#define G4VPGeoDriver_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PTouchableKey.hh"

class G4VPGeoDriver
{

public:  // with description

  virtual ~G4VPGeoDriver(){}
  
  virtual G4PTouchableKey LocateOnBoundary(const G4ThreeVector &aPosition, 
		                           const G4ThreeVector &aDirection) = 0;
    // Must be called in the PostDOIT of the ParallelTransportation.
    // The track crosses the boundary if PostDOIT gets called.
  
  virtual G4PTouchableKey GetCurrentTouchableKey() const = 0;

  virtual G4double ComputeStepLengthInit(const G4ThreeVector &aPosition, 
                                         const G4ThreeVector &aDirection) = 0;
  
  virtual G4double ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
				                  const G4ThreeVector &aDirection) = 0;

  virtual G4double ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
                                             const G4ThreeVector &aDirection) = 0;

};

#endif
