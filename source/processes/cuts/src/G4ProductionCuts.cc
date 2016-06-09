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
// $Id: G4ProductionCuts.cc,v 1.4 2004/06/07 13:47:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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


G4int  G4ProductionCuts::GetIndex(const G4String& name)
{
  G4int index;
  if       ( name == "gamma" )        { index =  0; }
  else  if ( name == "e-" )           { index =  1; }
  else  if ( name == "e+" )           { index =  2; }
  else                                { index = -1; }

  return index;
}


G4int  G4ProductionCuts::GetIndex(const G4ParticleDefinition* ptcl)
{ 
  if(!ptcl) return -1;
  if(gammaDef==0 && ptcl->GetParticleName()=="gamma") { gammaDef = ptcl; }
  if(electDef==0 && ptcl->GetParticleName()=="e-") { electDef = ptcl; }
  if(positDef==0 && ptcl->GetParticleName()=="e+") { positDef = ptcl; }
  G4int index;
  if(ptcl==gammaDef)      { index = 0; }
  else if(ptcl==electDef) { index = 1; }
  else if(ptcl==positDef) { index = 2; }
  else                    { index = -1; }

  return index;
}

