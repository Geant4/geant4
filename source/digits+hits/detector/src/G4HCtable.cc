// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCtable.cc,v 1.1 1999/01/07 16:06:26 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#include "G4HCtable.hh"

G4HCtable::G4HCtable() {;}

G4HCtable::~G4HCtable() {;}

G4int G4HCtable::Registor(G4String SDname,G4String HCname)
{
  for(int i=0;i<HClist.entries();i++)
  { if(HClist[i]==HCname && SDlist[i]==SDname) return -1; }
  HClist.insert(HCname);
  SDlist.insert(SDname);
  return HClist.entries();
}

G4int G4HCtable::GetCollectionID(G4String HCname)
{
  G4int i = -1;
  if(HCname.index("/")==RW_NPOS) // HCname only
  {
    for(G4int j=0;j<HClist.entries();j++)
    {
      if(HClist[j]==HCname)
      { 
        if(i>=0) return -2;
        i = j;
      }
    }
  }
  else
  {
    for(G4int j=0;j<HClist.entries();j++)
    {
      G4String tgt = SDlist[j];
      tgt += "/";
      tgt += HClist[j];
      if(tgt==HCname)
      {
        if(i>=0) return -2;
        i = j;
      }
    }
  }
  return i;
}


