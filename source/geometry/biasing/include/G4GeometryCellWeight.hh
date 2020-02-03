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
// G4GeometryCellWeight
//
// Class description:
//
// Used internally by weight window technique sampling. 
// It is a map of geometry cells to maps of upper energy
// to lower weight bounds.

// Author: Michael Dressel (CERN), 2002
// ----------------------------------------------------------------------
#ifndef G4GEOMETRYCELLWEIGHT_HH
#define G4GEOMETRYCELLWEIGHT_HH 1

#include <map>
#include "globals.hh"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

using G4UpperEnergyToLowerWeightMap =
  std::map<G4double, G4double, std::less<G4double> >;
using G4GeometryCellWeight =
  std::map<G4GeometryCell, G4UpperEnergyToLowerWeightMap, G4GeometryCellComp>;

#endif
