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
// $Id: G4MaterialCutsCouple.cc,v 1.2 2002-12-16 11:15:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    18 Sep. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4MaterialCutsCouple.hh"
#include "g4std/iomanip"

#include "G4Material.hh"
#include "G4ProductionCuts.hh"

G4MaterialCutsCouple::G4MaterialCutsCouple() :
  isMaterialModified(false),
  fMaterial(0),
  fCuts(0)
{
}
  
G4MaterialCutsCouple::G4MaterialCutsCouple(const G4Material* material,
					   G4ProductionCuts* cut) :
  isMaterialModified(true),
  fMaterial(material),
  fCuts(cut)
{
}


G4MaterialCutsCouple::G4MaterialCutsCouple(const G4MaterialCutsCouple& right) 
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

  return *this;
}

void G4MaterialCutsCouple::SetProductionCuts(G4ProductionCuts* aCut)
{ fCuts = aCut; }

G4ProductionCuts* G4MaterialCutsCouple::GetProductionCuts() const
{ return fCuts; }







