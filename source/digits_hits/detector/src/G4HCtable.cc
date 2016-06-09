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
// $Id: G4HCtable.cc,v 1.3 2004/03/15 19:16:07 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4HCtable.hh"
#include "G4VSensitiveDetector.hh"

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

G4int G4HCtable::GetCollectionID(G4String HCname) const
{
  G4int i = -1;
  if(HCname.index("/")==std::string::npos) // HCname only
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

G4int G4HCtable::GetCollectionID(G4VSensitiveDetector* aSD) const
{
  if(aSD->GetNumberOfCollections()<1)
  {
    G4cerr << "Sensitive detector <" << aSD->GetName()
           << "> does not have a registered hits collection."
           << G4endl;
    return -1;
  }
  if(aSD->GetNumberOfCollections()>1)
  {
    G4cerr << "Sensitive detector <" << aSD->GetName()
           << "> has more than one registered hits collections."
           << G4endl;
    G4cerr << "Candidates are : ";
    for(G4int j=0;j<aSD->GetNumberOfCollections();j++)
    { G4cerr << aSD->GetCollectionName(j) << " "; }
    G4cerr << G4endl;
    return -1;
  }
  for(size_t k=0;k<SDlist.size();k++)
  { if(SDlist[k]==aSD->GetName()) return k; }
  return -1;
}


