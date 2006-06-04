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
// $Id: pymodG4materials.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4materials.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4Version.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Material();
void export_G4MaterialTable();
void export_G4Element();
void export_G4ElementTable();
void export_G4NistManager();

BOOST_PYTHON_MODULE(G4materials)
{
  export_G4Material();
  export_G4MaterialTable();
  export_G4Element();
  export_G4ElementTable();
#if G4VERSION_NUMBER >= 710
  export_G4NistManager();
#endif
}

