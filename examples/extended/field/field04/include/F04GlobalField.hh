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

#ifndef F04GlobalField_h
#define F04GlobalField_h 1

#include <vector>

#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4MagneticField.hh"
#include "G4ElectroMagneticField.hh"

#include "G4Mag_EqRhs.hh"
#include "G4Mag_SpinEqRhs.hh"

#include "G4EqMagElectricField.hh"
#include "G4EqEMFieldWithSpin.hh"

#include "F04FieldMessenger.hh"
#include "F04ElementField.hh"

//  F04GlobalField - handles the global ElectroMagnetic field
//
//  There is a single G04GlobalField object.
//
//  The field from each individual beamline element is given by a
//  ElementField object. Any number of overlapping ElementField
//  objects can be added to the global field. Any element that
//  represents an element with an EM field must add the appropriate
//  ElementField to the global GlobalField object.

typedef std::vector<F04ElementField*> FieldList;

class F04GlobalField : public G4ElectroMagneticField {
//class F04GlobalField : public G4MagneticField {

private:

  F04GlobalField();
  F04GlobalField(const F04GlobalField&);

  ~F04GlobalField();

  F04GlobalField& operator=(const F04GlobalField&);

  void setupArray();

public:

  /// getObject() returns the single F04GlobalField object.
  /// It is constructed, if necessary.
  static F04GlobalField* getObject();

  /// GetFieldValue() returns the field value at a given point[].
  /// field is really field[6]: Bx,By,Bz,Ex,Ey,Ez.
  /// point[] is in global coordinates: x,y,z,t.
  void GetFieldValue(const G4double* point, G4double* field) const;

  /// DoesFieldChangeEnergy() returns true.
  G4bool DoesFieldChangeEnergy() const { return true; }

  /// addElementField() adds the ElementField object for a single
  /// element to the global field.
  void addElementField(F04ElementField* f) { if (fields) fields->push_back(f); }

  /// clear() removes all ElementField-s from the global object,
  /// and destroys them. Used before the geometry is completely
  /// re-created.
  void clear();

  /// updates all field tracking objects and clear()
  void updateField();

  /// Set the Stepper types
  void SetStepperType( G4int i ) { fStepperType = i; }

  /// Set the Stepper
  void SetStepper();

  /// Set the minimum step length
  void SetMinStep(G4double s) { minStep = s; }

  /// Set the delta chord length
  void SetDeltaChord(G4double s) { deltaChord = s; }

  /// Set the delta one step length
  void SetDeltaOneStep(G4double s) { deltaOneStep = s; }

  /// Set the delta intersection length
  void SetDeltaIntersection(G4double s) { deltaIntersection = s; }

  /// Set the minimum eps length
  void SetEpsMin(G4double s) { epsMin = s; }

  /// Set the maximum eps length
  void SetEpsMax(G4double s) { epsMax = s; }

  /// Return the list of Element Fields
  FieldList* getFields() { return fields; }

protected:

  /// Get the global field manager
  G4FieldManager* GetGlobalFieldManager();

private:

  static F04GlobalField* object;

  G4int nfp;
  G4bool first;

  FieldList* fields;

  const F04ElementField **fp;

private:

  G4int fStepperType;

  G4double minStep;
  G4double deltaChord;
  G4double deltaOneStep;
  G4double deltaIntersection;
  G4double epsMin;
  G4double epsMax;

//  G4Mag_EqRhs*            fEquation;
//  G4Mag_SpinEqRhs*        fEquation;

//  G4EqMagElectricField*   fEquation;
  G4EqEMFieldWithSpin*    fEquation;

  G4FieldManager*         fFieldManager;
  G4PropagatorInField*    fFieldPropagator;
  G4MagIntegratorStepper* fStepper;
  G4ChordFinder*          fChordFinder;

  F04FieldMessenger*         fFieldMessenger;

};

#endif
