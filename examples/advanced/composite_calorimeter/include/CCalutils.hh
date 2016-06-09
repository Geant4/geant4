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
///////////////////////////////////////////////////////////////////////////////
// File: CCalutils.hh
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalutils_hh
#define CCalutils_hh 1

#include <iostream>
#include <fstream>
#include "globals.hh"


G4String operator+(const G4String&, const int);
G4String operator+(const G4String&, const double);
// "number " + i = "number i"

std::ifstream& readName(std::ifstream&, G4String&);
// It reads a name into G4String between quotes and skips lines 
// beginning with #. and if found *ENDDO returns.

std::ifstream& findDO(std::ifstream&, const G4String&);
// It reads until a *DO str is found.

std::ostream& tab(std::ostream&);
// It add a tab.

std::istream& jump(std::istream&);
// It ignores character until the end of line.

bool openGeomFile(std::ifstream& is, const G4String& pathname, 
		  const G4String& filename);
// It opens the geometry file, either locally (if it exists) or "remotely".


#endif
