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
// $Id: pymodG4global.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4global.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_globals();
void export_geomdefs();
void export_G4StateManager();
void export_G4ApplicationState();
void export_G4String();
void export_G4ThreeVector();
void export_G4RotationMatrix();
void export_G4Transform3D();
void export_G4UnitsTable();
void export_Randomize();
void export_RandomEngines();
void export_G4RandomDirection();
void export_G4UserLimits();
void export_G4Timer();

BOOST_PYTHON_MODULE(G4global) 
{
  export_globals();
  export_geomdefs();
  export_G4StateManager();
  export_G4ApplicationState();
  export_G4String();
  export_G4ThreeVector();
  export_G4RotationMatrix();
  export_G4Transform3D();
  export_G4UnitsTable();
  export_Randomize();
  export_RandomEngines();
  export_G4RandomDirection();
  export_G4UserLimits();
  export_G4Timer();
}

