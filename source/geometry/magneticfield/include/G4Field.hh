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
// $Id: G4Field.hh,v 1.9 2004/09/01 09:59:06 japost Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
// class G4Field
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
// -------------------------------------------------------------------
// History:
// - Created:  John Apostolakis, 10.03.1997
// - Modified: 
//   V. Grichine     8 Nov 2001: Extended "Point" arg to [4] array to add time
//   J. Apostolakis  5 Nov 2003: Added virtual method DoesFieldChangeEnergy()
//   J. Apostolakis 31 Aug 2004: Information on convention for components
// -------------------------------------------------------------------

#ifndef G4FIELD_HH
#define G4FIELD_HH

#include "G4Types.hh"

class G4Field
{
  public:  // with description

      virtual void  GetFieldValue( const  double Point[4],
                                          double *fieldArr ) const = 0;
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

      G4Field(){;}
      virtual ~G4Field(){;}
      inline G4Field& operator = (const G4Field &p); 

     // A field signature function that can be used to insure
     // that the Equation of motion object and the G4Field object
     // have the same "field signature"?
 
      virtual G4bool   DoesFieldChangeEnergy() const= 0 ;
       //  Each type/class of field should respond this accordingly
       //  For example:
       //    - an electric field     should return "true"
       //    - a pure magnetic field should return "false"

};

inline G4Field& G4Field::operator = (const G4Field &p)
{
  if (&p != this) { *this = p; }
  return *this;
}

#endif /* G4FIELD_HH */
