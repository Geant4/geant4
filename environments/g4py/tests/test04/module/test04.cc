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
// $Id: test04.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   test04.cc
//
//   test of "Call Policies"
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "Step.hh"
#include "Track.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test04) {
  class_<Step>( "Step", "step class")
    .def(init<>())
    .add_property("x", &Step::GetX, &Step::SetX)
    ;

  class_<Track>( "Track", "track class")
    .def(init<>())
    .def("GetStep1", &Track::GetStep1,
	 return_internal_reference<>())
    .def("GetStep2", &Track::GetStep2,
	 return_value_policy<reference_existing_object>())
    // this is invalid, just for test
    .def("GetStep3", &Track::GetStep3,
	 return_value_policy<manage_new_object>())
    ;
}

