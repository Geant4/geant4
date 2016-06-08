// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCtable.cc,v 1.1.8.1.2.1 1999/12/07 20:47:44 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
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
  if(HCname.index("/")==G4std::string::npos) // HCname only
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


