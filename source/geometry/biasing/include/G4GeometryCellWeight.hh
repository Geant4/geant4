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
// $Id: G4GeometryCellWeight.hh,v 1.2 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// ----------------------------------------------------------------------
// Class G4GeometryCellWeight
//
// Class description:
//
// Used internally by weight window tehnique sampling. 
// It is a map of geometry cells to maps of upper energy
// to lower weight bounds.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCellWeight_hh
#define G4GeometryCellWeight_hh G4GeometryCellWeight_hh 

#include <map>
#include "globals.hh"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"
typedef std::map<G4double, G4double, std::less<G4double> > G4UpperEnergyToLowerWeightMap;
typedef std::map<G4GeometryCell, G4UpperEnergyToLowerWeightMap, G4GeometryCellComp>  G4GeometryCellWeight;

#endif
