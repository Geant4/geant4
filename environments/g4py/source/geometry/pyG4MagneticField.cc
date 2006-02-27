// $Id: pyG4MagneticField.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

