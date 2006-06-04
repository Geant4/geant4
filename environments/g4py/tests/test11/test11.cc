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
// $Id: test11.cc,v 1.3 2006-06-04 21:36:00 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test11.cc
//
//   test for the case for lack of default constructor
//
//                                         2005 Q
// ====================================================================
class AClass {
private:
  AClass() { }

public:
  AClass(int i)  { }
  AClass(const AClass& a) { }
  ~AClass() { }
};


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;


BOOST_PYTHON_MODULE(test11)
{
  class_<AClass>("AClass", "a class", no_init)
    .def(init<int>())
    ;
}

