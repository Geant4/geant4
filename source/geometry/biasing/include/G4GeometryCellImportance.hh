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
//
// $Id: G4GeometryCellImportance.hh,v 1.3 2003/06/16 16:50:59 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// Class G4GeometryCellImportance
//
// Class description:
//
// Used internally by importance sampling. 
// It is acontainer for "cell" importance value pairs.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCellImportance_hh
#define G4GeometryCellImportance_hh G4GeometryCellImportance_hh 

#include <map>
#include "globals.hh"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

typedef std::map<G4GeometryCell, G4double, G4GeometryCellComp> G4GeometryCellImportance;
  // implement container G4GeometryCellImportance as map

std::ostream& operator<<(std::ostream &out, const G4GeometryCellImportance &gCelli);

#endif
