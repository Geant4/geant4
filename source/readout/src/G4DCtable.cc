// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCtable.cc,v 1.3 1999-12-15 14:53:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DCtable.hh"

G4DCtable::G4DCtable() {;}

G4DCtable::~G4DCtable() {;}

G4int G4DCtable::Registor(G4String DMname,G4String DCname)
{
  for(int i=0;i<DClist.entries();i++)
  { if(DClist[i]==DCname && DMlist[i]==DMname) return -1; }
  DClist.insert(DCname);
  DMlist.insert(DMname);
  return DClist.entries();
}

G4int G4DCtable::GetCollectionID(G4String DCname)
{
  G4int i = -1;
  if(DCname.index("/")==G4std::string::npos) // DCname only
  {
    for(G4int j=0;j<DClist.entries();j++)
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
    for(G4int j=0;j<DClist.entries();j++)
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


