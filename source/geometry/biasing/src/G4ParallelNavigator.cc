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
// $Id: G4ParallelNavigator.cc,v 1.7 2002-07-29 15:56:20 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelNavigator.cc
//
// ----------------------------------------------------------------------

#include "g4std/strstream"

#include "G4ParallelNavigator.hh"

#include "G4Navigator.hh"
#include "G4PTouchableKey.hh"
#include "G4VParallelStepper.hh"

G4ParallelNavigator::G4ParallelNavigator(G4VPhysicalVolume &aWorldVolume)
  : fNavigator(*(new G4Navigator)), 
  fNlocated(0),
  fMaxShiftedTrys(10)  
{
  fNavigator.SetWorldVolume(&aWorldVolume);
  fCurrentTouchableH = fNavigator.CreateTouchableHistory();
}

G4ParallelNavigator::~G4ParallelNavigator()
{
  delete &fNavigator;
}

// public functions

G4PTouchableKey G4ParallelNavigator::
LocateOnBoundary(const G4ThreeVector &aPosition, 
		 const G4ThreeVector &aDirection)
{/*
  fNavigator.SetGeometricallyLimitedStep();
  fNavigator.
    LocateGlobalPointAndUpdateTouchableHandle( aPosition,
					       aDirection,
					       fCurrentTouchableH,
					       true);
  fNlocated++;
 */
  fNavigator.SetGeometricallyLimitedStep();
  Locate(aPosition, aDirection, false, true);
  //  since the track crosses a boundary ipdate stepper 
  return GetCurrentTouchableKey();
}

G4double G4ParallelNavigator::
ComputeStepLengthInit(const G4ThreeVector &aPosition, 
		      const G4ThreeVector &aDirection)
{
  // initialization
  fNlocated = 0;

  Locate(aPosition, aDirection, false, false);
  return GetStepLength("ComputeStepLengthInit",
			 aPosition, aDirection);
}

G4double G4ParallelNavigator::
ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
			       const G4ThreeVector &aDirection)
{
  if (fNlocated == 0 ) {
    Error("ComputeStepLengthCrossBoundary: no location done before call",
	  aPosition, aDirection);
  }

  // if the track is on boundary the LocateOnBoundary was called
  // in the DOIT of the ParalleTransport


  return  GetStepLength("ComputeStepLengthCrossBoundary",
			 aPosition, aDirection);
}


G4double G4ParallelNavigator::
ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
                          const G4ThreeVector &aDirection)
{
  if (fNlocated == 0 ) {
    Error("ComputeStepLengthInVolume, no location done before call", 
	  aPosition, aDirection);
  }

  // if the track is not on the boundary and it's
  // not the first step, the location must be inside the 
  // volume 

  fNavigator.LocateGlobalPointWithinVolume(aPosition);
  //  Locate(aPosition, aDirection, true); // reduces stepLength<=0. calls
  return  GetStepLengthUseLocate("ComputeStepLengthInVolume",
				 aPosition, aDirection);
}

// private functions

void G4ParallelNavigator::Locate(const G4ThreeVector &aPosition, 
                                 const G4ThreeVector &aDirection,
				 G4bool historysearch,
				 G4bool useDirection)
{

  fNavigator.LocateGlobalPointAndSetup( aPosition, 
					&aDirection, 
					historysearch, 
					!useDirection);
  fCurrentTouchableH = fNavigator.CreateTouchableHistory();

  fNlocated++;
  
  return;
}


G4double G4ParallelNavigator::
GetStepLength(const G4String methodname,
	      const G4ThreeVector &aPosition, 
	      const G4ThreeVector &aDirection) {
  G4double newSafety;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
						kInfinity, newSafety);
  // if stepLength = 0 try shifting
  if (stepLength<=0.) {
    stepLength  = 
      ComputeStepLengthShifted(methodname,
			       aPosition, aDirection);
  }
  return stepLength;
}

G4double G4ParallelNavigator::
GetStepLengthUseLocate(const G4String methodname,
	      const G4ThreeVector &aPosition, 
	      const G4ThreeVector &aDirection) {
  G4double newSafety;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
						kInfinity, newSafety);
  // if stepLength = 0 try locate
  if (stepLength<=0.) {
    Locate(aPosition, aDirection, true, true);    
    stepLength  = GetStepLength(methodname + 
				"form GetStepLengthUseLocate",
				aPosition, 
				aDirection);
  }
  return stepLength;
}


G4double
G4ParallelNavigator::
ComputeStepLengthShifted(const G4String &m,
			 const G4ThreeVector &aPosition, 
			 const G4ThreeVector &aDirection)
{
  G4cerr << "G4ParallelNavigator::ComputeStepLengthShifted: invoked by: "
	 << m << G4endl;
  G4ThreeVector shift_pos(aPosition);
  shift_pos+=G4ThreeVector(Shift(aDirection.x()), 
			   Shift(aDirection.y()), 
			   Shift(aDirection.z()));
  G4double stepLength = 0.;
  G4int trys = 0;
  while (stepLength<=0 && trys < fMaxShiftedTrys) {
    G4cerr << "G4ParallelNavigator::ComputeStepLengthShifted: trys = "
	   << ++trys << G4endl;
    G4cerr << "shifted position: " << shift_pos 
	   << ", direction: " << aDirection << G4endl;
    Locate(shift_pos, aDirection, true, true);
    G4double newSafety;
    stepLength = fNavigator.ComputeStep( shift_pos, aDirection,
					 kInfinity, newSafety);
  }
  if (stepLength<=0.0) {
    G4std::ostrstream os;
    os << "still got stepLength<=0.0: " << shift_pos << '\0' << G4endl;
    Error(os.str() , aPosition, aDirection);
  } 
  return stepLength;
}

G4PTouchableKey G4ParallelNavigator::GetCurrentTouchableKey() const
{
  return G4PTouchableKey(*fCurrentTouchableH->GetVolume(),
			 fCurrentTouchableH->GetReplicaNumber());
}

void G4ParallelNavigator::Error(const G4String &m,
                                const G4ThreeVector &pos,
                                const G4ThreeVector &dir)
{
  G4cerr << "ERROR - G4ParallelNavigator::" << m << G4endl;
  G4cerr << "aPosition: " << pos << G4endl;
  G4cerr << "dir: " << dir << G4endl;
  G4Exception("Program aborted.");
}

G4double G4ParallelNavigator::Shift(G4double d)
{
  if (d>0) return 2 * kCarTolerance;
  if (d<0) return -2 * kCarTolerance;
  return 0.;
}
