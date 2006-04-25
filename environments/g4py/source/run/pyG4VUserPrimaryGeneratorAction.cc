// $Id: pyG4VUserPrimaryGeneratorAction.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VUserPrimaryGeneratorAction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VUserPrimaryGeneratorAction {

struct CB_G4VUserPrimaryGeneratorAction :
  G4VUserPrimaryGeneratorAction, wrapper<G4VUserPrimaryGeneratorAction> {
  
  void GeneratePrimaries(G4Event* anEvent) {
    get_override("GeneratePrimaries")(&anEvent);
  }
};

};

using namespace pyG4VUserPrimaryGeneratorAction;

// ====================================================================
// module definition
// ====================================================================
void export_G4VUserPrimaryGeneratorAction()
{
  class_<CB_G4VUserPrimaryGeneratorAction, boost::noncopyable> 
    ("G4VUserPrimaryGeneratorAction", 
     "base class of user primary generator action")

    .def("GeneratePrimaries", 
	 pure_virtual(&G4VUserPrimaryGeneratorAction::GeneratePrimaries))
    ;
}

