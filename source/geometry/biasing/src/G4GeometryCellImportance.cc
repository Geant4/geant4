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
//
// $Id: G4GeometryCellImportance.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryCellImportance.cc
//
// ----------------------------------------------------------------------

#include "G4GeometryCellImportance.hh"
#include "G4GeometryCellStepStream.hh"

std::ostream& operator<<(std::ostream &out, const G4GeometryCellImportance &gCelli)
{
  for (G4GeometryCellImportance::const_iterator it = gCelli.begin();
       it != gCelli.end(); ++it) {
    out << (*it).first << ", importance = ";
    out << (*it).second << "\n";
  }
  return out;
}
