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
//   G4MCTSimParticle.cc
//
// ====================================================================
#include "g4std/strstream"
#include "g4std/iomanip"
#include "G4ios.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4MCTSimParticle.hh"
#include "G4MCTSimVertex.hh"

// ====================================================================
//
// class description
//
// ====================================================================

/////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle()
  : parentParticle(0), 
    trackID(0), parentTrackID(0),
    primaryFlag(false), 
    vertex(0), storeFlag(false)
/////////////////////////////////
{
}

/////////////////////////////////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle(G4std::string aname, int apcode, 
			       int atid, int ptid,
			       const HepLorentzVector& p)
  : parentParticle(0), 
    name(aname), pdgID(apcode), 
    trackID(atid), parentTrackID(ptid),
    primaryFlag(false),momentumAtVertex(p),
    vertex(0), storeFlag(false)
///////////////////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle(G4std::string aname, int apcode, 
			       int atid, int ptid,
			       const HepLorentzVector& p,
			       const G4MCTSimVertex* v )
  : parentParticle(0), 
    name(aname), pdgID(apcode), 
    trackID(atid), parentTrackID(ptid),
    primaryFlag(false),momentumAtVertex(p), 
    vertex(const_cast<G4MCTSimVertex*>(v)), storeFlag(false)
///////////////////////////////////////////////////////////////
{
}

/////////////////////////////////
G4MCTSimParticle::~G4MCTSimParticle()
/////////////////////////////////
{
  associatedParticleList.clear();
}

////////////////////////////////////////////////////////
int G4MCTSimParticle::AssociateParticle(G4MCTSimParticle* p)
////////////////////////////////////////////////////////
{
  associatedParticleList.push_back(p);
  p-> SetParentParticle(this);
  return associatedParticleList.size();
}

/////////////////////////////////////////////////////
int G4MCTSimParticle::GetNofAssociatedParticles() const
/////////////////////////////////////////////////////
{
  return associatedParticleList.size();
}

//////////////////////////////////////////////////////////////////
G4MCTSimParticle* G4MCTSimParticle::GetAssociatedParticle(int i) const
//////////////////////////////////////////////////////////////////
{
  int size= associatedParticleList.size();
  if(i>=0 && i< size) return associatedParticleList[i];
  else return 0;
}

////////////////////////////////////////
int G4MCTSimParticle::GetTreeLevel() const
////////////////////////////////////////
{
  const G4MCTSimParticle* p= this;
  int nlevel;
  for(nlevel=1;;nlevel++) {
    p= p-> GetParentParticle();
    if(p==0) return nlevel;
  }
}

/////////////////////////////////////////////////////
void G4MCTSimParticle::SetStoreFlagToParentTree(G4bool q)
/////////////////////////////////////////////////////
{
  storeFlag=q;
  if(vertex) vertex-> SetStoreFlag(q);
  if(primaryFlag) return;
  if(parentParticle) parentParticle-> SetStoreFlagToParentTree(q);
}


//////////////////////////////////////////////////////////
void G4MCTSimParticle::PrintSingle(G4std::ostream& ostr) const
//////////////////////////////////////////////////////////
{
  char stp[20];
  G4std::ostrstream os(stp,20);
  char cqp=' ';
  if(storeFlag) cqp='+';
  os << cqp << trackID << '\0';
  G4std::string stid(stp);
  ostr << G4std::setw(6) << stid;
  //ostr << G4std::setw(4) << trackID;

  if(primaryFlag) ostr << "*";
  else ostr << " ";
  ostr << "<" << G4std::setw(5) << parentTrackID;
  ostr.setf(G4std::ios::fixed);
  ostr << ": P(" 
      << G4std::setw(7) << G4std::setprecision(3) << momentumAtVertex.x()/GeV 
      << "," << G4std::setw(7) << G4std::setprecision(3) 
      << momentumAtVertex.y()/GeV  
      << "," << G4std::setw(7) << G4std::setprecision(3) 
      << momentumAtVertex.z()/GeV 
      << "," << G4std::setw(7) << G4std::setprecision(3) 
      << momentumAtVertex.e()/GeV << ") @";
  ostr << name << "(" << pdgID << ")";

  if(vertex) {
    ostr << " %" << vertex-> GetCreatorProcessName() << G4endl;

    char stv[20];
    G4std::ostrstream osv(stv,20);
    char cqv=' ';
    if(vertex->GetStoreFlag()) cqv='+';
    osv << cqv << vertex-> GetID() << '\0';
    G4std::string svid(stv);
    ostr << "       " << G4std::setw(6) << svid;
    //ostr << "      " << G4std::setw(4) << vertex-> GetID();
    ostr.unsetf(G4std::ios::fixed);
    ostr.setf(G4std::ios::scientific|G4std::ios::right|G4std::ios::showpoint);
    ostr << "- X(" << G4std::setw(9) << G4std::setprecision(2) 
	<< vertex-> GetPosition().x()/mm 
	<< "," << G4std::setw(9) << G4std::setprecision(2) 
	<< vertex-> GetPosition().y()/mm
	<< "," << G4std::setw(9) << G4std::setprecision(2) 
	<< vertex-> GetPosition().z()/mm 
	<< "," << G4std::setw(9) << G4std::setprecision(2) 
	<< vertex-> GetTime()/ns << ")";
    ostr.unsetf(G4std::ios::scientific);
    
    ostr << " @" << vertex-> GetVolumeName()
	<< "-" << vertex-> GetVolumeNumber();
  } 
  ostr << G4endl;
  
}

////////////////////////////////////////////////////////////////////
void G4MCTSimParticle::Print(G4std::ostream& ostr, G4bool qrevorder) const
////////////////////////////////////////////////////////////////////
{
  PrintSingle(ostr);

  // recursively print associated particles
  if (!qrevorder) { // parent -> child
    SimParticleList::const_iterator itr;
    for(itr= associatedParticleList.begin(); 
	itr!= associatedParticleList.end(); ++itr) {
      (*itr)-> Print(ostr);
    }
  } else { // child -> parent
    if(parentParticle) parentParticle-> Print(ostr, true);
  }
}
