///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Ecal.hh
// Description: Equipped to describe crystal matrix for different testbeam run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Ecal_h
#define CCalG4Ecal_h 1

#include "CCalEcal.hh"
#include "CCalG4Able.hh"
#include <vector>

typedef G4LogicalVolume* ptrG4Log;

class CCalG4Ecal: public CCalEcal, public CCalG4Able {
public:
  //Backward or Forward type
  enum CMType {module1, module2};

  //Constructor and Destructor
  CCalG4Ecal(const G4String &name);
  virtual ~CCalG4Ecal();

  void setType(CMType ty)    {type = ty;}
  
  //Prefix to all names in the Detector
  static G4String idName;  

protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() ;

private:
  //Methods to construct the different parts of the detector
  G4LogicalVolume* constructGlobal();

private:
  //Private data members
  CMType type;

  //Static logical volumes shared by forward and backward detectors.
  static G4LogicalVolume* crystalmatrixLog;

  // Logical volumes for sensitive detectors
  vector<ptrG4Log> sensitiveLogs;

};

#endif














