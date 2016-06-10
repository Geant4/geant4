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
// $Id: pyG4MagneticField.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4MagneticField.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================

class PyG4MagneticField : public G4MagneticField {
public:
  PyG4MagneticField() { }
  ~PyG4MagneticField() { }

  virtual G4ThreeVector GetFieldValue(const G4ThreeVector& pos,
				      const G4double time) const = 0;

  virtual void GetFieldValue(const G4double Point[4],
			     G4double* Bfield) const {

    const G4ThreeVector& bfield=
      GetFieldValue(G4ThreeVector(Point[0], Point[1], Point[2]), Point[3]);

    Bfield[0]= bfield.x();
    Bfield[1]= bfield.y();
    Bfield[2]= bfield.z();
  }

};

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4MagneticField {

struct CB_PyG4MagneticField :
    PyG4MagneticField, wrapper<PyG4MagneticField> {

  G4ThreeVector GetFieldValue(const G4ThreeVector& pos,
			      const G4double time) const {
    return get_override("GetFieldValue")(pos, time);
  }

};

G4ThreeVector(PyG4MagneticField::*f1_GetFieldValue)
  (const G4ThreeVector&, const G4double) const
  = &PyG4MagneticField::GetFieldValue;

}

using namespace pyG4MagneticField;

// ====================================================================
// module definition
// ====================================================================
void export_G4MagneticField()
{
  class_<G4MagneticField, boost::noncopyable >
    ("__G4MagneticField", "dummy class of magnetic field", no_init)
    ;

  class_<CB_PyG4MagneticField, boost::noncopyable,
    bases<G4Field, G4MagneticField> >
    ("G4MagneticField", "base class of magnetic field")
    // ---
    .def("DoesFieldChangeEnergy", &G4MagneticField::DoesFieldChangeEnergy)
    .def("GetFieldValue", pure_virtual(f1_GetFieldValue))
    ;
}

