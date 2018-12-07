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
/// \file field/field04/include/F04FocusSolenoid.hh
/// \brief Definition of the F04FocusSolenoid class
//

#ifndef F04FocusSolenoid_h
#define F04FocusSolenoid_h 1

#include "G4LogicalVolume.hh"

#include "F04SimpleSolenoid.hh"

//  F04FocusSolenoid implements a solenoid magnet with an increased 
//  cylindrical magnetic field in +/- z-direction (or half focusing)

class F04FocusSolenoid : public F04SimpleSolenoid
{
  public:

    ///  Default constructor.
    F04FocusSolenoid(G4double, G4double, G4double,
                            G4LogicalVolume*, G4ThreeVector);

    ///  Destructor.
    virtual ~F04FocusSolenoid() {}

    ///  Set F04FocusSolenoid to only half-focusing (default negative z)
    void SetHalf(G4bool h) { fHalf = h; }

    ///  AddFieldValue() adds the field for this solenoid into field[].
    ///  point[] is in global coordinates.
    virtual void AddFieldValue(const G4double point[4], G4double field[6]) const;

  private:

  G4bool fHalf;

  G4double fB1;
  G4double fB2;

};

#endif
