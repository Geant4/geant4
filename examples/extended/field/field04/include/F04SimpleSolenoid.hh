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
//
/// \file field/field04/include/F04SimpleSolenoid.hh
/// \brief Definition of the F04SimpleSolenoid class
//

#ifndef F04SimpleSolenoid_h
#define F04SimpleSolenoid_h 1

#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"

#include "F04ElementField.hh"
#include "F04GlobalField.hh"

//  F04SimpleSolenoid implements a solenoid magnet with a constant
//  cylindrical magnetic field in z-direction

class F04SimpleSolenoid : public F04ElementField
{
  public:

    ///  Default constructor.
    F04SimpleSolenoid(G4double, G4double, G4LogicalVolume*, G4ThreeVector);

    ///  Destructor.
    virtual ~F04SimpleSolenoid() {}

    ///  GetLength() returns the length of the solenoid
    virtual G4double GetLength() { return fFieldLength; }

    ///  GetWidth() returns the solenoid diameter
    virtual G4double GetWidth() { return fFieldRadius*2.0; }

    ///  GetHeight() returns the solenoid diameter
    virtual G4double GetHeight() { return fFieldRadius*2.0; }

    /// SetFringeZ(G4double) sets the solenoid fringe field z-length
    void SetFringeZ(G4double z) { fFringeZ = z; }

    /// GetFringeZ() returns the solenoid fringe field z-length
    G4double GetFringeZ() { return fFringeZ; }

    ///  IsOutside() returns true when outside the solenoid
    G4bool IsOutside(G4ThreeVector& local) const;

    /// IsWithin() returns true when inside the solenoid
    G4bool IsWithin(G4ThreeVector& local) const;

    ///  AddFieldValue() adds the field for this solenoid into field[].
    ///  point[] is in global coordinates.
    virtual void AddFieldValue(const G4double point[4], G4double field[6]) const;

  private:

    G4double fBfield;

    G4double fFringeZ;

    G4double fFieldLength;
    G4double fFieldRadius;

};

#endif
