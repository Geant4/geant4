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
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

    ///  getLength() returns the length of the solenoid
    virtual G4double getLength() { return fieldLength; }

    ///  getWidth() returns the solenoid diameter
    virtual G4double getWidth() { return fieldRadius*2.0; }

    ///  getHeight() returns the solenoid diameter
    virtual G4double getHeight() { return fieldRadius*2.0; }

    /// setFringeZ(G4double) sets the solenoid fringe field z-length
    void setFringeZ(G4double z) { fringeZ = z; }

    /// getFringeZ() returns the solenoid fringe field z-length
    G4double getFringeZ() { return fringeZ; }

    ///  isOutside() returns true when outside the solenoid
    G4bool isOutside(G4ThreeVector& local) const; 

    /// isWithin() returns true when inside the solenoid
    G4bool isWithin(G4ThreeVector& local) const;

    ///  addFieldValue() adds the field for this solenoid into field[].
    ///  point[] is in global coordinates.
    void addFieldValue(const G4double point[4], G4double field[6]) const;

  private:

    G4double Bfield;

    G4double fringeZ;

    G4double fieldRadius;
    G4double fieldLength;

};

#endif
