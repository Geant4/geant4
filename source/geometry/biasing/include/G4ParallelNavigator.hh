#ifndef G4ParallelNavigator_hh
#define G4ParallelNavigator_hh G4ParallelNavigator_hh

#include "G4VPGeoDriver.hh"
#include "geomdefs.hh"


class G4VTouchable;
class G4Navigator;
class G4PTouchableKey ;
class G4VPhysicalVolume;

class G4ParallelNavigator : public G4VPGeoDriver {
public:
  G4ParallelNavigator(G4VPhysicalVolume &aWorldVolume);
  ~G4ParallelNavigator();
  
  // from G4VPGeoDriver
  G4PTouchableKey 
  LocateOnBoundary(const G4ThreeVector &aPosition, 
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
			
  void Error(const G4String &m, const G4ThreeVector &pos, const G4ThreeVector &dir);
  

  G4double Shift(G4double d) {
    if (d>0) return 2 * kCarTolerance;
    if (d<0) return -2 * kCarTolerance;
    return 0;
  }
  
  G4Navigator &fNavigator;
  G4VTouchable *fCurrentTouchable;
  G4int fNlocated;

};

#endif













