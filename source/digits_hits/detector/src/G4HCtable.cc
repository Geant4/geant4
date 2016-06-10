//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4HCtable.cc 67992 2013-03-13 10:59:57Z gcosmo $
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


