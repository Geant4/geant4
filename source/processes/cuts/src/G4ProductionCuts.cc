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
// $Id: G4ProductionCuts.cc,v 1.2 2003/09/19 14:45:27 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    18 Sep. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4ProductionCuts.hh"
#include <iomanip>

const G4ParticleDefinition* G4ProductionCuts::gammaDef = 0;
const G4ParticleDefinition* G4ProductionCuts::electDef = 0;
const G4ParticleDefinition* G4ProductionCuts::positDef = 0;

G4ProductionCuts::G4ProductionCuts() :
  isModified(true)
{
  for (G4int i=0; i< NumberOfG4CutIndex; i++) {
    fRangeCuts.push_back(0.0);
  }
}

G4ProductionCuts::G4ProductionCuts(const G4ProductionCuts& right) 
{
  *this = right;
}

G4ProductionCuts::~G4ProductionCuts()
{
  fRangeCuts.clear();
}

G4ProductionCuts & G4ProductionCuts::operator=(const G4ProductionCuts &right)
{
  if (&right==this) return *this;

  for (G4int i=0; i< NumberOfG4CutIndex; i++) {
    fRangeCuts[i] = right.fRangeCuts[i];
  }
  isModified = right.isModified;

  return *this;
}



G4int G4ProductionCuts::operator==(const G4ProductionCuts &right) const
{
  return (this == &right);
}


G4int G4ProductionCuts::operator!=(const G4ProductionCuts &right) const
{
  return (this !=  &right);
}

