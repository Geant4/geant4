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
// $Id: G4ParallelNavigator.cc,v 1.6 2002-07-18 14:55:50 dressel Exp $
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
  : fNavigator(*(new G4Navigator)), fNlocated(0)
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
{
  fNavigator.SetGeometricallyLimitedStep();
  Locate(aPosition, aDirection, true);
  // since the track crosses a boundary ipdate stepper 
  return GetCurrentTouchableKey();
}

G4double G4ParallelNavigator::
ComputeStepLengthInit(const G4ThreeVector &aPosition, 
		      const G4ThreeVector &aDirection)
{
  G4double newSafety;
  G4double stepLength;

  // initialization
  fNlocated = 0;
  
  // first try
  Locate(aPosition, aDirection, true);
  stepLength  = fNavigator.ComputeStep( aPosition, aDirection,
					kInfinity, newSafety);
  if (stepLength<=0.0) {
    // second try with schifted position
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthInit",
					   aPosition, aDirection);
  } 
  return stepLength;
}

G4double G4ParallelNavigator::
ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
			       const G4ThreeVector &aDirection)
{
  if (fNlocated == 0 ) {
    Error("ComputeStepLengthCrossBoundary: no location done before call",
	  aPosition, aDirection);
  }

  // if the track is on boundary the location was done
  // in the DOIT of the ImportanceSplitExaminer

  G4double newSafety;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
				       kInfinity, newSafety);
  if (stepLength<=0.) {
    // second try with schifted position
    G4cout << "Warning: G4ParallelNavigator::ComputeStepLengthCrossBoundary: "
	   << G4endl;
    G4cout << "stepLength<=0." << ", pos = " << aPosition 
	   << ", dir = " << aDirection << G4endl;
    G4cout << "try ComputeStepLengthShifted" << G4endl;
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthCrossBoundary",
					   aPosition, aDirection);
  }
  return stepLength;
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

  //  fNavigator.LocateGlobalPointWithinVolume(aPosition);
  Locate(aPosition, aDirection, true); // reduces stepLength<=0. calls

  G4double newSafety;
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
				       kInfinity, newSafety);
  if (stepLength<=0.) {
    // second try with schifted position
    G4cout << "Warning: G4ParallelNavigator::ComputeStepLengthInVolume: "
	   << G4endl;
    G4cout << "stepLength<=0." << ", pos = " << aPosition 
	   << ", dir = " << aDirection << G4endl;
    G4cout << "try ComputeStepLengthShifted" << G4endl;
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthInVolume",
					   aPosition, aDirection);
  }
  return stepLength;
}

// private functions

void G4ParallelNavigator::Locate(const G4ThreeVector &aPosition, 
                                 const G4ThreeVector &aDirection,
                                       G4bool useDirection)
{

  fNavigator.LocateGlobalPointAndSetup( aPosition, &aDirection, true, !useDirection);
  fCurrentTouchableH = fNavigator.CreateTouchableHistory();

  fNlocated++;
  
  return;
}

G4double
G4ParallelNavigator::ComputeStepLengthShifted(const G4String &m,
                                              const G4ThreeVector &aPosition, 
                                              const G4ThreeVector &aDirection)
{
  G4ThreeVector shift_pos(aPosition);
  shift_pos+=G4ThreeVector(Shift(aDirection.x()), 
			   Shift(aDirection.y()), 
			   Shift(aDirection.z()));
  Locate(shift_pos, aDirection, true);
  G4double newSafety;
  G4double stepLength = fNavigator.ComputeStep( shift_pos, aDirection,
                                                kInfinity, newSafety);
  if (stepLength<=0.0) {
    G4std::ostrstream os;
    os << "ComputeStepLengthShifted: called by " << m << "\n";
    os << "shifted positio in second try: " << shift_pos << '\0' << G4endl;
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
  G4cout << "ERROR - G4ParallelNavigator::" << m << G4endl;
  G4cout << "aPosition: " << pos << G4endl;
  G4cout << "dir: " << dir << G4endl;
  G4Exception("Program aborted.");
}

G4double G4ParallelNavigator::Shift(G4double d)
{
  if (d>0) return 2 * kCarTolerance;
  if (d<0) return -2 * kCarTolerance;
  return 0.;
}
