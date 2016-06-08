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
