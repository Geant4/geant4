#ifndef G4VPGeoDriver_hh 
#define  G4VPGeoDriver_hh G4VPGeoDriver_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PTouchableKey.hh"



class G4VPGeoDriver {
public:
  virtual ~G4VPGeoDriver(){}
  
  // LocateOnBoundary must be called in the PostDOIT
  // of the ParallelTransportation. The track crosses
  // the boundary if PostDOIT get's called.
  virtual G4PTouchableKey 
  LocateOnBoundary(const G4ThreeVector &aPosition, 
		   const G4ThreeVector &aDirection) = 0;  
  
  
  virtual G4PTouchableKey GetCurrentTouchableKey() const = 0;

  virtual G4double
  ComputeStepLengthInit(const G4ThreeVector &aPosition, 
			const G4ThreeVector &aDirection) = 0;
  
  virtual G4double 
  ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
				 const G4ThreeVector &aDirection) = 0;
  
  virtual G4double 
  ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
			    const G4ThreeVector &aDirection) = 0;

};

#endif
