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
// $Id: test06.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   test06.cc
//
//   test for inheritance
//   - inherit C++ base class in Python side
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "XBase.hh"
#include "ZBase.hh"

using namespace boost::python;

struct w_XBase : XBase, wrapper<XBase> {
  std::string PVMethod() {
    return this-> get_override("PVMethod")();
  }
};


struct w_ZBase : ZBase, wrapper<ZBase> {
  void VMethod(std::string message) {
    if(override f= this-> get_override("VMethod"))
      f(message);
    else 
      ZBase::VMethod(message);
  }

  void d_VMethod(std::string message) {
    this-> ZBase::VMethod(message);
  }
};


BOOST_PYTHON_MODULE(test06) {
  //class_<w_XBase, boost::noncopyable>("XBase", "base class", no_init)
  class_<w_XBase, boost::noncopyable>("XBase", "base class")
    .add_property("ival", &XBase::GetIVal, &XBase::SetIVal)
    .def("PVMethod", pure_virtual(&XBase::PVMethod))
    ;

  class_<w_ZBase, boost::noncopyable>("ZBase", "base class")
    .def("AMethod", &ZBase::AMethod)
    .def("VMethod", &ZBase::VMethod, &w_ZBase::d_VMethod)
    ;
}
