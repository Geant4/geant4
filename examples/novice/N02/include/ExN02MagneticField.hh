// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02MagneticField.hh,v 1.3 2000-12-04 16:24:05 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#ifndef ExN02MagneticField_H
#define ExN02MagneticField_H

#include "G4UniformMagField.hh"

class G4FieldManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


class ExN02MagneticField: public G4UniformMagField
{
  public:
  
   ExN02MagneticField(G4ThreeVector);  //  The value of the field
   ExN02MagneticField();               //  A zero field
  ~ExN02MagneticField();  
      
   //Set the field (fieldValue,0,0)
   void SetFieldValue(G4double fieldValue);
   void SetFieldValue(G4ThreeVector fieldVector);
      
   G4ThreeVector GetConstantFieldValue();

  protected:

      // Find the global Field Manager
      G4FieldManager* GetGlobalFieldManager();   // static 
};

#endif
