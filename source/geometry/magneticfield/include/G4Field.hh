//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4Field
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.
//
// The key method is  GetFieldValue( const  double Point[4],
//                    *************         double *fieldArr ) 
// Given an input position/time vector 'Point', 
// this method must return the value of the field in "fieldArr".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy() 
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, optionally, 
// spin.  For this a field and its equation of motion must follow the
// same convention for the order of field components in the array "fieldArr"

// Created: John Apostolakis, 10.03.1997
// -------------------------------------------------------------------
#ifndef G4FIELD_HH
#define G4FIELD_HH

#include "G4Types.hh"
#include "globals.hh"

class G4Field
{
  public:  // with description

      G4Field( G4bool gravityOn = false);
      G4Field( const G4Field& );
      virtual ~G4Field();
      G4Field& operator = (const G4Field& p); 

      virtual void GetFieldValue( const G4double Point[4],
                                        G4double* fieldArr ) const = 0;
       // Given the position time vector 'Point', 
       // return the value of the field in the array fieldArr.
       //  Notes: 
       //   1) The 'Point' vector has the following structure:
       //        Point[0]  is  x  ( position, in Geant4 units )
       //        Point[1]  is  y
       //        Point[2]  is  z
       //        Point[3]  is  t  ( time,  in Geant4 units )
       //   2) The convention for the components of the field
       //      array 'fieldArr' are determined by the type of field.
       //      See for example the class G4ElectroMagneticField.

      virtual G4bool DoesFieldChangeEnergy() const = 0;
        // Each type/class of field should respond this accordingly
        // For example:
        //   - an electric field     should return "true"
        //   - a pure magnetic field should return "false"

      inline G4bool IsGravityActive() const;
        //  Does this field include gravity?

      inline void SetGravityActive( G4bool OnOffFlag );
    
      virtual G4Field* Clone() const;
        // Implements cloning, needed by multi-threading

      static constexpr G4int MAX_NUMBER_OF_COMPONENTS = 24;

  private:

      G4bool fGravityActive = false;
};

// Inline methods ...

inline G4bool G4Field::IsGravityActive() const
{
  return fGravityActive;
}

inline void G4Field::SetGravityActive( G4bool OnOffFlag )
{ 
  fGravityActive = OnOffFlag; 
}
 
#endif
