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
// $Id: G4ParallelNavigator.hh,v 1.2 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelNavigator
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelNavigator_hh
#define G4ParallelNavigator_hh

#include "G4VPGeoDriver.hh"
#include "geomdefs.hh"

class G4VTouchable;
class G4Navigator;
class G4PTouchableKey ;
class G4VPhysicalVolume;

class G4ParallelNavigator : public G4VPGeoDriver
{

public:  // with description

  G4ParallelNavigator(G4VPhysicalVolume &aWorldVolume);
  ~G4ParallelNavigator();
  
  G4PTouchableKey LocateOnBoundary(const G4ThreeVector &aPosition, 
		                   const G4ThreeVector &aDirection);

  G4PTouchableKey GetCurrentTouchableKey() const;

  G4double ComputeStepLengthInit(const G4ThreeVector &aPosition, 
				 const G4ThreeVector &aDirection);
  
  G4double ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
					  const G4ThreeVector &aDirection);
  
  G4double ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
				     const G4ThreeVector &aDirection);
  
private:

  G4ParallelNavigator(const G4ParallelNavigator &);
  G4ParallelNavigator &operator=(const G4ParallelNavigator &);

  G4double ComputeStepLengthShifted(const G4String &m,
				    const G4ThreeVector &aPosition, 
				    const G4ThreeVector &aDirection);

  void Locate(const G4ThreeVector &aPosition, 
	      const G4ThreeVector &aDirection,
	      G4bool histsearch); 
			
  void Error(const G4String &m,
             const G4ThreeVector &pos,
             const G4ThreeVector &dir);

  G4double Shift(G4double d);

private:
  
  G4Navigator &fNavigator;
  G4VTouchable *fCurrentTouchable;
  G4int fNlocated;
};

#endif













