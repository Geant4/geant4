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
// $Id: G4DCtable.cc,v 1.6 2001-10-18 20:28:19 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DCtable.hh"

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

G4int G4DCtable::GetCollectionID(G4String DCname)
{
  G4int i = -1;
  if(DCname.index("/")==G4std::string::npos) // DCname only
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


