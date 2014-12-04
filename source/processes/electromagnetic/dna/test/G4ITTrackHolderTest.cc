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
/*
 * G4ITTrackHolderTest.cc
 *
 *  Created on: 17 nov. 2014
 *      Author: kara
 */

#include "G4ITTrackHolder.hh"
#include "G4Electron_aq.hh"
#include "G4OH.hh"
#include "G4H2O2.hh"
#include "G4Hydrogen.hh"
#include "G4Molecule.hh"
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4MoleculeTable.hh"
#include "G4MoleculeFinder.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"

using namespace std;

void Iterate()
{
  G4TrackManyList* manyLists = G4ITTrackHolder::Instance()->GetMainList();
  G4TrackManyList::iterator it = manyLists->begin();
  G4TrackManyList::iterator end = manyLists->end();

  G4cout << "Iterate" << G4endl;

  for (; it != end; it++)
  {
    if (*it == 0)
    {
      G4cerr << "No track found " << G4endl;
      abort();
    }

    G4cout << "trackID : " << (*it)->GetTrackID();
    G4cout << " " << GetMolecule(*it)->GetName()
    << G4endl;
  }

}

void Clear(vector<G4Molecule*> vec)
{
  for (int i = 0; i < vec.size(); i++)
  {
    if (vec[i])
    {
      delete vec[i]->GetTrack();
      vec[i] = 0;
    }
  }

}

void PushTrack(G4Track* track)
{
  G4ITTrackHolder::Instance()->Push(track);
  G4MoleculeFinder::Instance()->Push(track);
}


vector<G4ThreeVector> OH_position;
vector<G4ThreeVector> OHm_position;


double FindClosest(const G4ThreeVector& _from, vector<G4ThreeVector>& _in)
{
  size_t _out_i = 0;
  double dis = 3e8;

  for(size_t i = 0 ; i != _in.size() ; i++)
  {
    double tmp= (_in[i]-_from).mag();
    if(tmp < dis)
    {
      _out_i = i;
      dis = tmp;
    }
  }

  return dis;
}

//pair<G4ThreeVector,G4ThreeVector>
double
FindClosest(vector<G4ThreeVector>& _from, vector<G4ThreeVector>& _in)
{
  size_t _out_i = 0;
  double dis = 3e8;

  for(size_t i = 0 ; i != _from.size() ; i++)
  {
    double tmp= FindClosest(_from[i],_in);
    if(tmp < dis)
    {
      _out_i = i;
      dis = tmp;
    }
  }

  return dis;
}

double
FindClosest(vector<G4ThreeVector>& _in)
{
  size_t _out_i = 0;
  double dis = 3e8;

  for(size_t i = 0 ; i != _in.size() ; i++)
  {
    for(size_t j = i+1 ; j != _in.size() ; j++)
    {
      double tmp= (_in[i]-_in[j]).mag();
      if(tmp < dis)
      {
        _out_i = i;
        dis = tmp;
      }
    }
  }

  return dis;
}


int main()
{
  G4MoleculeTable::Instance()->CreateMoleculeModel("H2O2",
                                                   G4H2O2::Definition());
  G4MoleculeTable::Instance()->CreateMoleculeModel("H",
                                                   G4Hydrogen::Definition());
  G4MoleculeTable::Instance()->CreateMoleculeModel("OHm",
                                                   G4OH::Definition(),1);

  std::vector<G4Molecule*> e_aq;
  std::vector<G4Molecule*> OH;
  std::vector<G4Molecule*> OHm;
  std::vector<G4Molecule*> H2O2;
  std::vector<G4Molecule*> H;

  for (int i = 0; i < 10; i++)
  {
    e_aq.push_back(new G4Molecule(G4Electron_aq::Definition()));
    G4Track* track =  e_aq[i]->BuildTrack(1 * picosecond, G4ThreeVector());
    PushTrack(track);
  }

  for (int i = 0; i < 10; i++)
  {
    OH.push_back(new G4Molecule(G4OH::Definition()));
    G4ThreeVector pos = G4RandomDirection()*G4UniformRand()*100;
    G4Track* track =  OH[i]->BuildTrack(1 * picosecond, pos);
    OH_position.push_back(pos);
    PushTrack(track);
  }

  double time;
  G4ITTrackHolder::Instance()->MergeNextTimeToMainList(time);

  G4MoleculeFinder::Instance()->UpdatePositionMap();
  G4KDTreeResultHandle kdtree_res = G4MoleculeFinder::Instance()->FindNearest(OH[0],OH[0]);
  double kdtree_dis = sqrt(kdtree_res->GetDistanceSqr());

  G4cout << "kdtree_dis" << kdtree_dis << G4endl;
  G4cout << "expected dis" << FindClosest(OH_position) << G4endl;

  assert(fabs(kdtree_dis - FindClosest(OH_position)) < 1e-6);

  for (int i = 0; i < 1; i++)
  {
    H2O2.push_back(new G4Molecule(G4H2O2::Definition()));
    G4Track* track =  H2O2[i]->BuildTrack(1 * picosecond, G4ThreeVector());
    PushTrack(track);
  }

  for (int i = 0; i < 1; i++)
  {
    H.push_back(new G4Molecule(G4Hydrogen::Definition()));
    G4Track* track =  H[i]->BuildTrack(1 * picosecond, G4ThreeVector());
    PushTrack(track);
  }

  G4ITTrackHolder::Instance()->MergeNextTimeToMainList(time);

  G4MoleculeFinder::Instance()->UpdatePositionMap();
  kdtree_res = G4MoleculeFinder::Instance()->FindNearest(OH[0],OH[0]);
  kdtree_dis = sqrt(kdtree_res->GetDistanceSqr());

  assert(fabs(kdtree_dis - FindClosest(OH_position)) < 1e-6);

  std::map<Key, PriorityList*>& priorityLists = G4ITTrackHolder::Instance()
      ->GetLists();
  assert(priorityLists.size() == 4);

  G4cout << "-------------------" << G4endl;
  G4cout << "Test priority list" << G4endl;
  {
    std::map<Key, PriorityList*>::iterator it = priorityLists.begin();
    std::map<Key, PriorityList*>::iterator end = priorityLists.end();
    for (; it != end; it++)
    {
      assert(it->second->Get(PriorityList::MainList));
      assert(it->second->Get(PriorityList::MainList)->size() != 0);

      assert(it->second->Get(PriorityList::WaitingList) == 0);
      //assert(it->second->Get(PriorityList::WaitingList)->size() == 0);

      assert(it->second->Get(PriorityList::SecondariesList));
      assert(it->second->Get(PriorityList::SecondariesList)->size() == 0);
    }

  }
  G4cout << "Priority list OK" << G4endl;

  G4cout << "-------------------" << G4endl;
  G4cout << "Test number of tracks" << G4endl;
  G4cout << "-------------------" << G4endl;
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 22);
  G4cout << "Number of tracks OK" << G4endl;

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << "Test first iteration" << G4endl;
  G4cout << "-------------------" << G4endl;
  {
    std::map<Key, PriorityList*>::iterator it = priorityLists.begin();
    std::map<Key, PriorityList*>::iterator end = priorityLists.end();

    for (; it != end; it++)
    {
      G4TrackList* trackList = it->second->Get(PriorityList::MainList);

      if(trackList->empty())
      {
        G4cout << "Is empty" << G4endl;
      }

      G4TrackList::iterator it2 = trackList->begin();
      G4TrackList::iterator end2 = trackList->end();
      for (; it2 != end2; it2++)
      {
        G4cout << GetMolecule(*it2)->GetName() << G4endl;
      }
    }
  }
  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << "Test second iteration" << G4endl;
  G4cout << "-------------------" << G4endl;
  Iterate();
  G4cout << "Iteration OK" << G4endl;

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << "Withdraw one e_aq" << G4endl;
  G4cout << "-------------------" << G4endl;
  delete e_aq[5]->GetTrack();
  e_aq[5] = 0;
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 21);
  Iterate();
  G4cout << "Withdraw one e_aq OK" << G4endl;

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Clear one H2O2" << G4endl;
  G4cout << "-------------------" << G4endl;
  delete H2O2[0]->GetTrack();
  H2O2[0] = 0;
  PriorityList* H2O2_list = priorityLists[G4MoleculeTable::Instance()
      ->GetMoleculeModel("H2O2")->GetMoleculeID()];

  assert(H2O2_list);
  assert(H2O2_list->GetMainList()->size() == 0);

  assert(G4ITTrackHolder::Instance()->GetNTracks() == 20);

  assert(priorityLists.begin()->second->GetMainList()->empty() == true);
  assert(priorityLists.begin()->second->GetNTracks() == 0);
  Iterate();

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Clear all H" << G4endl;
  G4cout << "-------------------" << G4endl;
  Clear(H);
  G4cout << G4ITTrackHolder::Instance()->GetNTracks() << G4endl;
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 19);

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Test adding a H" << G4endl;
  G4cout << "-------------------" << G4endl;
  PriorityList* H_list = priorityLists[G4MoleculeTable::Instance()
      ->GetMoleculeModel("H")->GetMoleculeID()];
  assert(H_list->GetMainList()->size() == 0);

  H[0] = new G4Molecule(G4Hydrogen::Definition());
  G4Track* track = H[0]->BuildTrack(1 * picosecond, G4ThreeVector());
  PushTrack(track);
  G4ITTrackHolder::Instance()->MergeNextTimeToMainList(time);

  assert(H_list->GetMainList()->size() == 1);
  assert(H_list->Get(PriorityList::WaitingList) == 0);
  assert(H_list->Get(PriorityList::SecondariesList)->size() == 0);

  bool found = false;
  G4TrackManyList* manyLists = G4ITTrackHolder::Instance()->GetMainList();
  G4TrackManyList::iterator it = manyLists->begin();
  G4TrackManyList::iterator end = manyLists->end();

  G4cout << "Iterate" << G4endl;

  size_t i = 0;
  for (; it != end; it++, i++)
  {
    G4cout << "i: " << i;
    G4cout << " "
    << GetMolecule(*it)->GetName() << " "<< (*it)->GetTrackID() << G4endl;
    if (GetMolecule(*it) == H[0])
    {
      found = true;
      break;
    }
  }

  assert(found);
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 20);

  G4cout << "-------------------" << G4endl;
  G4cout << " Test removing the unique H" << G4endl;
  G4cout << "-------------------" << G4endl;

  Clear(H);
  Iterate();

  G4cout << G4ITTrackHolder::Instance()->GetNTracks() << G4endl;
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 19);

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Clear all e_aq" << G4endl;
  G4cout << "-------------------" << G4endl;
  Clear(e_aq);
  Iterate();

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Clear OH" << G4endl;
  G4cout << "-------------------" << G4endl;
  Clear(OH);
  Iterate();

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << " Clear H2O2" << G4endl;
  G4cout << "-------------------" << G4endl;
  Clear(H2O2);
  Iterate();

  G4cout << G4endl;
  G4cout << "-------------------" << G4endl;
  G4cout << "Test final number of tracks = "
  << G4ITTrackHolder::Instance()->GetNTracks() << G4endl;
  assert(G4ITTrackHolder::Instance()->GetNTracks() == 0);
  G4cout << "-------------------" << G4endl;
  G4ITTrackHolder::Instance()->Clear();
  G4cout << "Test OK" << G4endl;
  return 0;
}

