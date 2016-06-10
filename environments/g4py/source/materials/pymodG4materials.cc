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
// $Id: pymodG4materials.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pymodG4materials.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Material();
void export_G4MaterialTable();
void export_G4Element();
void export_G4ElementTable();
//void export_G4ElementVector();
void export_G4Isotope();
void export_G4NistManager();
void export_G4AtomicShells();

BOOST_PYTHON_MODULE(G4materials)
{
  export_G4Material();
  export_G4MaterialTable();
  export_G4Element();
  export_G4ElementTable();
  //export_G4ElementVector();
  export_G4Isotope();
  export_G4NistManager();
  export_G4AtomicShells();
}

