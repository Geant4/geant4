// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08MagneticField.hh,v 1.1 1999-01-08 16:35:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
// 
//  $ Id:  $

// Should this be a:
//    i) messenger
//   ii) user class that creates the field       ? 
//  iii) simply a derived class of Uniform field ?  <== I have chosen this now.
//   iv) a field manager that creates/updates field    (Prefered?)
#ifndef T08MagneticField_H
#define T08MagneticField_H

#include "G4UniformMagField.hh"
class G4FieldManager;

class T08MagneticField: public G4UniformMagField
{
  public:
      T08MagneticField(G4ThreeVector);  //  The value of the field
      T08MagneticField();               //  A zero field
      ~T08MagneticField();  
      
      // Set the field to (0, 0, fieldValue)
      void SetFieldValue(G4ThreeVector fieldVector);
      void SetFieldValue(G4double      fieldValue);
      G4ThreeVector GetConstantFieldValue();

  protected:

      // Find the global Field Manager
      G4FieldManager* GetGlobalFieldManager();   // static 
};

#endif
