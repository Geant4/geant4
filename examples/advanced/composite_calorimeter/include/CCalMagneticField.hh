///////////////////////////////////////////////////////////////////////////////
// File: CCalMagneticField.hh
// Description: A class for control of the Magnetic Field of the detector.
//              The field is assumed to be uniform.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalMagneticField_H
#define CCalMagneticField_H

#include "G4UniformMagField.hh"
#include "G4ThreeVector.hh"
class G4FieldManager;

class CCalMagneticField: public G4MagneticField {
public:
  CCalMagneticField(const G4String &name);
  ~CCalMagneticField();  
      
  // Access functions
  void MagneticField(const double Point[3], double Bfield[3]) const;
  Hep3Vector MagneticField(const Hep3Vector Point) const;
  virtual void GetFieldValue(const double Point[3], double* Bfield) const;
  G4double GetConstantFieldvalue() const {return fval;}

protected:
  // Find the global Field Manager
  G4FieldManager* GetGlobalFieldManager(); 

private:
  G4double   fval;             // Field value
  G4int      npts;             // Number of poinst
  G4double   xoff;             // Offset
  G4double*  pos;              // Position
  G4double*  slope;            // Slope
  G4double*  intercept;        // Intercept
};

#endif
