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
// $Id: G4HCtable.cc,v 1.6 2001/07/13 15:00:08 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-00 $
//

#include "G4HCtable.hh"

G4HCtable::G4HCtable() {;}

G4HCtable::~G4HCtable() {;}

G4int G4HCtable::Registor(G4String SDname,G4String HCname)
{
  for(size_t i=0;i<HClist.size();i++)
  { if(HClist[i]==HCname && SDlist[i]==SDname) return -1; }
  HClist.push_back(HCname);
  SDlist.push_back(SDname);
  return HClist.size();
}

G4int G4HCtable::GetCollectionID(G4String HCname)
{
  G4int i = -1;
  if(HCname.index("/")==G4std::string::npos) // HCname only
  {
    for(size_t j=0;j<HClist.size();j++)
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
    for(size_t j=0;j<HClist.size();j++)
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


