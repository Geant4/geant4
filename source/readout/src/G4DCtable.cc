// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCtable.cc,v 1.4 2001-02-08 06:07:21 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DCtable.hh"

G4DCtable::G4DCtable() {;}

G4DCtable::~G4DCtable() {;}

G4int G4DCtable::Registor(G4String DMname,G4String DCname)
{
  for(int i=0;i<DClist.size();i++)
  { if(DClist[i]==DCname && DMlist[i]==DMname) return -1; }
  DClist.push_back(DCname);
  DMlist.push_back(DMname);
  return DClist.size();
}

G4int G4DCtable::GetCollectionID(G4String DCname)
{
  G4int i = -1;
  if(DCname.index("/")==G4std::string::npos) // DCname only
  {
    for(G4int j=0;j<DClist.size();j++)
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
    for(G4int j=0;j<DClist.size();j++)
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


