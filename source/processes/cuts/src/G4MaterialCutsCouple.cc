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
// $Id: G4MaterialCutsCouple.cc 70369 2013-05-29 14:59:24Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    18 Sep. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4MaterialCutsCouple.hh"
#include <iomanip>

#include "G4Material.hh"
#include "G4ProductionCuts.hh"

G4MaterialCutsCouple::G4MaterialCutsCouple() :
  isMaterialModified(false),
  fMaterial(0),
  fCuts(0),
  indexNumber(-1),
  isUsedInGeometry(false)
{
}
  
G4MaterialCutsCouple::G4MaterialCutsCouple(const G4Material* material,
					   G4ProductionCuts* cut) :
  isMaterialModified(true),
  fMaterial(material),
  fCuts(cut),
  indexNumber(-1),
  isUsedInGeometry(false)
{
}


G4MaterialCutsCouple::G4MaterialCutsCouple(const G4MaterialCutsCouple& right) 
  :fMaterial(0), fCuts(0)
{
  *this = right;
}

G4MaterialCutsCouple::~G4MaterialCutsCouple()
{
}

G4MaterialCutsCouple & G4MaterialCutsCouple::operator=(const G4MaterialCutsCouple &right)
{
  if (&right==this) return *this;

  fMaterial = right.fMaterial;
  fCuts     = right.fCuts; 
  isMaterialModified = right.isMaterialModified;
  indexNumber = right.indexNumber;
  isUsedInGeometry = right.isUsedInGeometry;

  return *this;
}

