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
// $Id: test04.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
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

