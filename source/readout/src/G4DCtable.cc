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
// $Id: G4DCtable.cc 80987 2014-05-19 10:50:22Z gcosmo $
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

