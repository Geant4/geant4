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
// $Id: pyG4ProcVector.cc,v 1.4 2006-06-29 15:34:52 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ProcVector.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4indexing.hh"
#include "G4VProcess.hh"

using namespace boost::python;

namespace pyG4ProcVector {

typedef std::vector<G4VProcess*> G4ProcVector;

};

using namespace pyG4ProcVector;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcVector()
{
  class_<G4ProcVector> ("G4ProcVector", "process vector")
    .def(vector_indexing_suite<G4ProcVector>())
    ;
}

