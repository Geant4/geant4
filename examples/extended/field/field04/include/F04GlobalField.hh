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
/// \file field/field04/include/F04GlobalField.hh
/// \brief Definition of the F04GlobalField class
//

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

#include "F04DetectorConstruction.hh"

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

  F04GlobalField(F04DetectorConstruction* const);
  F04GlobalField(const F04GlobalField&);

  F04GlobalField& operator=(const F04GlobalField&);

  void SetupArray();

public:

  virtual ~F04GlobalField();

  /// GetObject() returns the single F04GlobalField object.
  /// It is constructed, if necessary.
  static F04GlobalField* GetObject(F04DetectorConstruction* const);
  static F04GlobalField* GetObject();

  /// GetFieldValue() returns the field value at a given point[].
  /// field is really field[6]: Bx,By,Bz,Ex,Ey,Ez.
  /// point[] is in global coordinates: x,y,z,t.
  virtual void GetFieldValue(const G4double* point, G4double* field) const;

  /// DoesFieldChangeEnergy() returns true.
  virtual G4bool DoesFieldChangeEnergy() const { return true; }

  /// AddElementField() adds the ElementField object for a single
  /// element to the global field.
  void AddElementField(F04ElementField* f)
  {
    if (fFields) fFields->push_back(f);
  }

  /// Clear() removes all ElementField-s from the global object,
  /// and destroys them. Used before the geometry is completely
  /// re-created.
  void Clear();

  /// constructs all field tracking objects
  void ConstructField();

  /// Set the Stepper types
  void SetStepperType( G4int i ) { fStepperType = i; }

  /// Set the Stepper
  void SetStepper();

  /// Set the minimum step length
  void SetMinStep(G4double stp) { fMinStep = stp; }

  /// Set the delta chord length
  void SetDeltaChord(G4double dcr) { fDeltaChord = dcr; }

  /// Set the delta one step length
  void SetDeltaOneStep(G4double stp) { fDeltaOneStep = stp; }

  /// Set the delta intersection length
  void SetDeltaIntersection(G4double its) { fDeltaIntersection = its; }

  /// Set the minimum eps length
  void SetEpsMin(G4double eps) { fEpsMin = eps; }

  /// Set the maximum eps length
  void SetEpsMax(G4double eps) { fEpsMax = eps; }

  /// Return the list of Element Fields
  FieldList* GetFields() { return fFields; }

protected:

  /// Get the global field manager
  G4FieldManager* GetGlobalFieldManager();

private:

  static G4ThreadLocal F04GlobalField* fObject;

  G4int fNfp;
  G4bool fFirst;

  FieldList* fFields;

  const F04ElementField **fFp;

private:

  G4int fStepperType;

  G4double fMinStep;
  G4double fDeltaChord;
  G4double fDeltaOneStep;
  G4double fDeltaIntersection;
  G4double fEpsMin;
  G4double fEpsMax;

//  G4Mag_EqRhs*            fEquation;
//  G4Mag_SpinEqRhs*        fEquation;

//  G4EqMagElectricField*   fEquation;
  G4EqEMFieldWithSpin*    fEquation;

  G4FieldManager*         fFieldManager;
  G4PropagatorInField*    fFieldPropagator;
  G4MagIntegratorStepper* fStepper;
  G4ChordFinder*          fChordFinder;

  F04FieldMessenger*         fFieldMessenger;

  F04DetectorConstruction* fDetectorConstruction;

};

#endif
