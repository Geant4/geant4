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
// $Id: G4ProductionCuts.cc,v 1.5 2003-04-02 21:56:47 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    18 Sep. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4ProductionCuts.hh"
#include "g4std/iomanip"

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
//  else  if ( name == "proton" )       { index =  3; }
//  else  if ( name == "anti_proton" )  { index =  4; }
//  else  if ( name == "neutron" )      { index =  5; }
//  else  if ( name == "anti_neutron" ) { index =  6; }
  else                                { index = -1; }

  return index;
}

G4int  G4ProductionCuts::GetIndex(const G4ParticleDefinition* ptcl)
{ return GetIndex(ptcl->GetParticleName()); }
  
void  G4ProductionCuts::SetProductionCut(G4double cut, G4int index)
{
  if (index<0) {
    for(G4int i = 0; i < NumberOfG4CutIndex; i++) {
      fRangeCuts[i] = cut;
    }
    isModified = true;

  } else if (index < NumberOfG4CutIndex) {
    fRangeCuts[index] = cut;
    isModified = true;
  }     
}

  
void  G4ProductionCuts::SetProductionCut(G4double cut, G4ParticleDefinition* ptcl)
{
  G4int idx = -1;
  if(ptcl) idx = GetIndex(ptcl->GetParticleName());
  if(idx>=0) SetProductionCut(cut,idx);
}


G4double  G4ProductionCuts::GetProductionCut(G4int index) const
{
  G4double cut=-1.0;
  if ( (index>=0) && (index<NumberOfG4CutIndex) ) {
    cut = fRangeCuts[index]; 
  }
  return cut;
}


G4double  G4ProductionCuts::GetProductionCut(const G4String& name) const
{
  return GetProductionCut(GetIndex(name)); 
}


void  G4ProductionCuts::SetProductionCuts(G4std::vector<G4double>& cut)
{
  for(G4int i = 0; (i<NumberOfG4CutIndex); i++) {
    fRangeCuts[i] = cut[i];
  }
  isModified = true;
}


const G4std::vector<G4double>&   G4ProductionCuts::GetProductionCuts() const
{
  return fRangeCuts;
}



G4bool  G4ProductionCuts::IsModified() const
{
  return isModified;
}


void   G4ProductionCuts::PhysicsTableUpdated()
{
  isModified = false;
}








