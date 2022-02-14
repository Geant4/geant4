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
// Authors: M. Large, S. Guatelli -  University of Wollongong, Australia
// This class developed by John Allison, March 2021
//

#ifndef ICRP110PhantomVisAction_h
#define ICRP110PhantomVisAction_h 1

#include "G4VUserVisAction.hh"

#include "G4Colour.hh"

class G4Material;

#include <set>

class ICRP110PhantomConstruction;

class ICRP110PhantomVisAction: public G4VUserVisAction
{
public:
  ICRP110PhantomVisAction(ICRP110PhantomConstruction* dc)
  : fpDetectorConstruction(dc) {}
  void Draw();
private:
  ICRP110PhantomConstruction*                    fpDetectorConstruction;
  std::multimap<const G4Material*,G4ThreeVector> fPositionByMaterial;
  std::map     <const G4Material*,G4Colour>      fColourByMaterial;
  std::set     <const G4Material*>               fSetOfMaterials;
};

#endif
