///////////////////////////////////////////////////////////////////////////////
// File: G4CrystalMatrix.hh
// Date:          06/09/9
// Modifications: 27/03/00 S.B. In OSCAR
// Description: Equipped to describe crystal matrix for different testbeam run
///////////////////////////////////////////////////////////////////////////////
#ifndef G4CrystalMatrix_h
#define G4CrystalMatrix_h 1

#include "CrystalMatrix.hh"
#include "CCalG4Able.hh"
#include <vector>

typedef G4LogicalVolume* ptrG4Log;

class G4CrystalMatrix: public CrystalMatrix, public CCalG4Able {
public:
  //Backward or Forward type
  enum CMType {module1, module2};

  //Constructor and Destructor
  G4CrystalMatrix(const G4String &name);
  virtual ~G4CrystalMatrix();

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














