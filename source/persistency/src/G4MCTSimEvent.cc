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
//   G4MCTSimEvent.cc
//
// ====================================================================
#include "G4ios.hh"
#include "G4MCTSimEvent.hh"
#include "G4MCTSimParticle.hh"
#include "G4MCTSimVertex.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////
G4MCTSimEvent::G4MCTSimEvent()
//////////////////////////
{
}

///////////////////////////
G4MCTSimEvent::~G4MCTSimEvent()
///////////////////////////
{
  ClearEvent();
}

//////////////////////////////////////////////////////////////
G4bool G4MCTSimEvent::AddParticle(const G4MCTSimParticle* aparticle)
//////////////////////////////////////////////////////////////
{
  G4MCTSimParticle* qpart= const_cast<G4MCTSimParticle*>(aparticle);
  int trackID= aparticle-> GetTrackID();
  int nc= particleMap.count(trackID);
  if(nc==0) {
    particleMap.insert(std::make_pair(trackID, qpart));
    return true;
  } else {
    return false;
  }

}

////////////////////////////////////////////////////////
G4MCTSimParticle* G4MCTSimEvent::FindParticle(int tid) const
////////////////////////////////////////////////////////
{
  G4MCTSimParticleContainer::const_iterator pos= particleMap.find(tid);
  if(pos != particleMap.end()) {
    return pos-> second;
  } else {
    return 0;
  }
}

///////////////////////////////////////////////////
G4MCTSimVertex* G4MCTSimEvent::GetVertex(int vid) const
///////////////////////////////////////////////////
{
  int nv= vertexVec.size();
  if(vid>=1 && vid<=nv) {
    return vertexVec[vid-1];
  } else {
    return 0;
  }
}

////////////////////////////////////////
void G4MCTSimEvent::BuildVertexContainer()
////////////////////////////////////////
{
  G4MCTSimParticleContainer::iterator itr;
  int vid=1;
  for(itr= particleMap.begin(); itr!= particleMap.end(); ++itr) {
    G4MCTSimVertex* vertex= itr->second-> GetVertex();    
    if(vertex) {      
      if (vertex-> GetID()<0) { // ID not yet assigned
	vertex-> SetID(vid);
	vid++;
	if (vertex) vertexVec.push_back(vertex);
      }
    }
  }
}

//////////////////////////////
void G4MCTSimEvent::ClearEvent() 
//////////////////////////////
{
  G4MCTSimParticleContainer::iterator itr;
  for(itr= particleMap.begin(); itr!= particleMap.end(); ++itr) {
    delete itr->second;
  }
  particleMap.clear();

  G4MCTSimVertexContainer::iterator itrv;
  for(itrv= vertexVec.begin(); itrv!= vertexVec.end(); ++itrv) {
    delete (*itrv);
  }
  vertexVec.clear();
}


//////////////////////////////////////////////
int G4MCTSimEvent::GetNofStoredParticles() const
//////////////////////////////////////////////
{
  int n=0;
  G4MCTSimParticleContainer::const_iterator itr;
  for(itr= particleMap.begin(); itr!= particleMap.end(); ++itr) {
    if(itr-> second->  GetStoreFlag()) n++;
  }
  return n;
}

/////////////////////////////////////////////
int G4MCTSimEvent::GetNofStoredVertices() const
/////////////////////////////////////////////
{
 int n=0;
  G4MCTSimVertexContainer::const_iterator itr;
  for(itr= vertexVec.begin(); itr!= vertexVec.end(); ++itr) {
    if((*itr)->GetStoreFlag()) n++;
  }
  return n;
}


/////////////////////////////////////////////////
void G4MCTSimEvent::Print(std::ostream& ostr) const
/////////////////////////////////////////////////
{
  ostr << "____________________________________________________"
          "____________________________" << G4endl;
  ostr << "SimEvent:" << G4endl << G4endl;
  ostr << "Current Memory Usage: " 
       << particleMap.size() << " particles, "
       << vertexVec.size() <<  " vertices."
       << G4endl;				      
  ostr << "trk#<ptrk#: P(Px(GeV),     Py,     Pz,     E ) @PDG     %proc\n"
       << "      vtx#- X(    X(mm),        Y,        Z,    T(ns)) @vname-#" 
       << G4endl;
  ostr << "____________________________________________________"
          "____________________________" << G4endl;

  G4MCTSimParticleContainer::const_iterator itr;
  for(itr= particleMap.begin(); itr!= particleMap.end(); ++itr) {
    itr-> second-> PrintSingle(ostr);
  }
  ostr << "____________________________________________________"
          "____________________________" << G4endl;
}
