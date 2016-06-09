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
// $Id: G4ParallelNavigator.cc,v 1.18 2003/10/22 13:52:10 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelNavigator.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelNavigator.hh"

#include <strstream>
#include "G4VParallelStepper.hh"


G4ParallelNavigator::G4ParallelNavigator(G4VPhysicalVolume &aWorldVolume)
  : 
  fNlocated(0),
  fMaxShiftedTrys(10),
  fCurrentTouchableH((fNavigator.SetWorldVolume(&aWorldVolume),
		      fNavigator.CreateTouchableHistory())),
  fVerbose(0)
{}

G4ParallelNavigator::~G4ParallelNavigator()
{
}

// public functions

G4GeometryCell G4ParallelNavigator::
LocateOnBoundary(const G4ThreeVector &aPosition, 
		 const G4ThreeVector &aDirection)
{/*
  fNavigator->SetGeometricallyLimitedStep();
  fNavigator->
    LocateGlobalPointAndUpdateTouchableHandle( aPosition,
					       aDirection,
					       fCurrentTouchableH,
					       true);
  fNlocated++;
 */
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "LOB: " << aPosition << ", " << aDirection << G4endl;
  }
  fNavigator.SetGeometricallyLimitedStep();
  Locate(aPosition, aDirection, true, false);
  //  since the track crosses a boundary ipdate stepper 
  G4GeometryCell g= GetCurrentGeometryCell();
  return g;
}

G4double G4ParallelNavigator::
ComputeStepLengthInit(const G4ThreeVector &aPosition, 
		      const G4ThreeVector &aDirection)
{
  // initialization
  fNlocated = 0;
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "CSInit: " << aPosition << ", " << aDirection << G4endl;
  }
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
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "CSCB: " << aPosition << ", " << aDirection << G4endl;
  }
  // if the track is on boundary the LocateOnBoundary was called
  // in the DOIT of the ParalleTransport

  // but location in volume overcomes a problem with spheres
  //  fNavigator->LocateGlobalPointWithinVolume(aPosition);
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
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "CSinV: " << aPosition << ", " << aDirection << G4endl;
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
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "Locate: " << aPosition << ", " << aDirection << G4endl;
  }
  fNavigator.LocateGlobalPointAndSetup( aPosition, 
					&aDirection, 
					historysearch, 
					!useDirection);
  fCurrentTouchableH = fNavigator.CreateTouchableHistory();

  fNlocated++;
  
  return;
}


G4double G4ParallelNavigator::
GetStepLength(const G4String &methodname,
	      const G4ThreeVector &aPosition, 
	      const G4ThreeVector &aDirection) {
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "GetSL: " << aPosition << ", " << aDirection << G4endl;
  }
  G4double newSafety = 0.;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, 
						aDirection,
						kInfinity, 
						newSafety);
  // if stepLength = 0 try shifting
  if (stepLength<=2*kCarTolerance) {
    stepLength  = 
      ComputeStepLengthShifted(methodname,
			       aPosition, aDirection);
  }
  return stepLength;
}

G4double G4ParallelNavigator::
GetStepLengthUseLocate(const G4String &methodname,
	      const G4ThreeVector &aPosition, 
	      const G4ThreeVector &aDirection) {
  if (fVerbose>=2) {
    G4cout.precision(12);
    G4cout << "GetSLuseLocate: " << aPosition << ", " << aDirection << G4endl;
  }
  G4double newSafety = 0;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
						kInfinity, newSafety);
  // if stepLength = 0 try locate
  if (stepLength<=2*kCarTolerance) {
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
  G4ThreeVector shift_pos(aPosition);
  G4double stepLength = 0.;
  G4int trys = 0;
  while (stepLength<=2*kCarTolerance && trys < fMaxShiftedTrys) {
    shift_pos+=G4ThreeVector(Shift(aDirection.x()), 
			     Shift(aDirection.y()), 
			     Shift(aDirection.z()));
    if (fVerbose>=1 || trys >= fMaxShiftedTrys - 1) {
      G4cout << "G4ParallelNavigator::ComputeStepLengthShifted: invoked by: "
	     << m << G4endl;
      G4cout << "  trys = "
	     << ++trys << G4endl;
      G4cout.precision(12);
      G4cout << "  shifted position: " << shift_pos 
	     << ", direction: " << aDirection << G4endl;
    }
    fNavigator.SetGeometricallyLimitedStep();
    Locate(shift_pos, aDirection, true, false); // to get into the correct volume
    G4double newSafety = 0;
    fNavigator.LocateGlobalPointWithinVolume(aPosition);    // to place at the correct position
    stepLength = fNavigator.ComputeStep( aPosition, aDirection,
					 kInfinity, newSafety);
  }
  if (stepLength<=kCarTolerance) {
    char st[200];
    std::ostrstream os(st,200);
    os << "still got stepLength<=kCarTolerance: "
       << shift_pos
       << "\n"
       << '\0';
    G4String m(st);
    Error(m, aPosition, aDirection);
  } 
  return stepLength;
}

G4GeometryCell G4ParallelNavigator::GetCurrentGeometryCell() const
{
  return G4GeometryCell(*fCurrentTouchableH->GetVolume(),
			 fCurrentTouchableH->GetReplicaNumber());
}


void  G4ParallelNavigator::SetVerboseity(G4int v) {
  fVerbose = v;
}

void G4ParallelNavigator::Error(const G4String &m,
                                const G4ThreeVector &pos,
                                const G4ThreeVector &dir)
{
  G4cerr << "ERROR - G4ParallelNavigator: " << m << G4endl
         << "        aPosition: " << pos << G4endl
         << "        dir: " << dir << G4endl;
  G4Exception("G4ParallelNavigator::Error()",
              "FatalException", FatalException, m);
}

G4double G4ParallelNavigator::Shift(G4double d)
{
  G4double s(0);
  G4double ds(10 * kCarTolerance);
  if (d>0){
    s = ds;
  }
  else if (d<0) {
    s = -1 * ds;
  }
  return s;
}
