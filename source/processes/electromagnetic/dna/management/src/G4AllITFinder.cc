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
// $Id: G4AllITFinder.cc 80074 2014-04-01 13:35:04Z matkara $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITFinder.hh"

using namespace std;

G4ThreadLocal G4AllITFinder* G4AllITFinder::fpInstance = 0;

G4AllITFinder::G4AllITFinder()
{
  fVerbose = 0;
}

G4AllITFinder* G4AllITFinder::Instance()
{
  if (!fpInstance) fpInstance = new G4AllITFinder();
  return fpInstance;
}

void G4AllITFinder::DeleteInstance()
{
  if (fpInstance)
  {
    delete fpInstance;
    fpInstance = 0;
  }
}

G4AllITFinder::~G4AllITFinder()
{
  std::map<G4ITType, G4VITFinder*>::iterator it;
  std::map<G4ITType, G4VITFinder*>::iterator it_tmp;

  for (it = fITSubManager.begin(); it != fITSubManager.end();)
  {
    if (it->second) delete it->second;
    it_tmp = it;
    it++;
    fITSubManager.erase(it_tmp);
  }
  fpInstance = 0;
}

void G4AllITFinder::UpdatePositionMap()
{
  std::map<G4ITType, G4VITFinder*>::iterator it = fITSubManager.begin();

  for (; it != fITSubManager.end(); it++)
  {
    it->second->UpdatePositionMap();
  }
}

G4VITFinder* G4AllITFinder::GetInstance(G4ITType type)
{
  map<G4ITType, G4VITFinder*>::iterator it = fITSubManager.find(type);

  if (it == fITSubManager.end()) return 0;

  return it->second;
}

void G4AllITFinder::RegisterManager(G4VITFinder* manager)
{
  fITSubManager[manager->GetITType()] = manager;
}

void G4AllITFinder::Push(G4Track* track)
{
  fITSubManager[GetIT(track)->GetITType()]->Push(track);
}

