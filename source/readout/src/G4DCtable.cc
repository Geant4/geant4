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
// $Id: G4DCtable.cc,v 1.8 2004/03/15 19:18:53 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4DCtable.hh"
#include "G4VDigitizerModule.hh"

G4DCtable::G4DCtable() {;}

G4DCtable::~G4DCtable() {;}

G4int G4DCtable::Registor(G4String DMname,G4String DCname)
{
  for(int i=0;i<int(DClist.size());i++)
  { if(DClist[i]==DCname && DMlist[i]==DMname) return -1; }
  DClist.push_back(DCname);
  DMlist.push_back(DMname);
  return DClist.size();
}

G4int G4DCtable::GetCollectionID(G4String DCname) const
{
  G4int i = -1;
  if(DCname.index("/")==std::string::npos) // DCname only
  {
    for(int j=0;j<int(DClist.size());j++)
    {
      if(DClist[j]==DCname)
      { 
        if(i>=0) return -2;
        i = j;
      }
    }
  }
  else
  {
    for(int j=0;j<int(DClist.size());j++)
    {
      G4String tgt = DMlist[j];
      tgt += "/";
      tgt += DClist[j];
      if(tgt==DCname)
      {
        if(i>=0) return -2;
        i = j;
      }
    }
  }
  return i;
}

G4int G4DCtable::GetCollectionID(G4VDigitizerModule* aDM) const
{
  if(aDM->GetNumberOfCollections()<1)
  {
    G4cerr << "Digitizer Module <" << aDM->GetName()
           << "> does not have a registered digits collection."
           << G4endl;
    return -1;
  }
  if(aDM->GetNumberOfCollections()>1)
  {
    G4cerr << "Digitizer Module <" << aDM->GetName()
           << "> has more than one registered digits collections."
           << G4endl;
    G4cerr << "Candidates are : ";
    for(G4int j=0;j<aDM->GetNumberOfCollections();j++)
    { G4cerr << aDM->GetCollectionName(j) << " "; }
    G4cerr << G4endl;
    return -1;
  }
  for(size_t k=0;k<DMlist.size();k++)
  { if(DMlist[k]==aDM->GetName()) return k; }
  return -1;
}

