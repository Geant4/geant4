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
// G4ElectroMagneticField
//
// Class description:
//
// A full electromagnetic field, containing both electric and magnetic fields.
// It is an abstract class, and a derived type of this field must be
// created by the user to describe his/her field configuration.
//
// We have established a convention for the electromagnetic field components:
// In the GetValue() method, the return values of Bfield will have 
// the following meaning
//  - Components 0, 1 and 2 are the Magnetic Field (x, y, z respectively);
//  - Components 3, 4 and 5 are the Electric field (x, y, z respectively).
// 
// Note 1: one or the other field could optional, depending on the Equation
// Note 2: such a convention is required between any field and its 
//         corresponding equation of motion.

// Created: J.Apostolakis, 12.11.1998
// Modified: V.Grichine, 08.11.2001: Extended "Point" to add time
// -------------------------------------------------------------------
#ifndef G4ELECTROMAGNETIC_FIELD_HH
#define G4ELECTROMAGNETIC_FIELD_HH

#include "G4Field.hh"

class G4ElectroMagneticField : public G4Field
{
  public:

    G4ElectroMagneticField();
   ~G4ElectroMagneticField() override;

    G4ElectroMagneticField(const G4ElectroMagneticField& r);
    G4ElectroMagneticField& operator = (const G4ElectroMagneticField& p);
      // Copy constructor & assignment operators.

    void  GetFieldValue(const G4double Point[4],
                              G4double *Bfield ) const override = 0;
      // Return as Bfield[0], [1], [2] the magnetic field x, y & z components
      // and    as Bfield[3], [4], [5] the electric field x, y & z components

    G4bool DoesFieldChangeEnergy() const override = 0;
      // For field with an electric component this should be true
      // For pure magnetic field this should be false
      // Alternative: default safe implementation { return true; }
};

#endif
