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
// $Id: pyG4MagneticField.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
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

};

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

