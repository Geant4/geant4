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
// $Id: pyG4AtomicShells.cc 101514 2016-11-18 15:30:57Z gcosmo $
// ====================================================================
//   pyG4AtomicShells.cc
//
//                                         2008 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4AtomicShells.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4AtomicShells()
{
  class_<G4AtomicShells, boost::noncopyable>
    ("G4AtomicShells", "Atomic subshell binding energy table", no_init)

    .def("GetNumberOfShells",    &G4AtomicShells::GetNumberOfShells)
    .staticmethod("GetNumberOfShells")

    .def("GetNumberOfElectrons", &G4AtomicShells::GetNumberOfElectrons)
    .staticmethod("GetNumberOfElectrons")

    .def("GetBindingEnergy", &G4AtomicShells::GetBindingEnergy)
    .staticmethod("GetBindingEnergy")

    .def("GetTotalBindingEnergy", &G4AtomicShells::GetTotalBindingEnergy)
    ;
}
