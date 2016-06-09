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
//   G4MCTSimParticle.cc
//
// ====================================================================

#include <sstream>
#include <iomanip>

#include "G4MCTSimParticle.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4MCTSimVertex.hh"

// ====================================================================
//
// class description
//
// ====================================================================

/////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle()
  : parentParticle(0), pdgID(0),
    trackID(0), parentTrackID(0),
    primaryFlag(false), 
    vertex(0), storeFlag(false)
/////////////////////////////////
{
}

/////////////////////////////////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle(std::string aname, int apcode, 
			       int atid, int ptid,
			       const G4LorentzVector& p)
  : parentParticle(0), 
    name(aname), pdgID(apcode), 
    trackID(atid), parentTrackID(ptid),
    primaryFlag(false),momentumAtVertex(p),
    vertex(0), storeFlag(false)
///////////////////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////////////////
G4MCTSimParticle::G4MCTSimParticle(std::string aname, int apcode, 
			       int atid, int ptid,
			       const G4LorentzVector& p,
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
void G4MCTSimParticle::PrintSingle(std::ostream& ostr) const
//////////////////////////////////////////////////////////
{
  std::ostringstream os;
  char cqp=' ';
  if(storeFlag) cqp='+';
  os << cqp << trackID << '\0';
  std::string stid(os.str());
  ostr << std::setw(6) << stid;
  //ostr << std::setw(4) << trackID;

  if(primaryFlag) ostr << "*";
  else ostr << " ";
  ostr << "<" << std::setw(5) << parentTrackID;
  ostr.setf(std::ios::fixed);
  ostr << ": P(" 
      << std::setw(7) << std::setprecision(3) << momentumAtVertex.x()/GeV 
      << "," << std::setw(7) << std::setprecision(3) 
      << momentumAtVertex.y()/GeV  
      << "," << std::setw(7) << std::setprecision(3) 
      << momentumAtVertex.z()/GeV 
      << "," << std::setw(7) << std::setprecision(3) 
      << momentumAtVertex.e()/GeV << ") @";
  ostr << name << "(" << pdgID << ")";

  if(vertex) {
    ostr << " %" << vertex-> GetCreatorProcessName() << G4endl;

    std::ostringstream osv;
    char cqv=' ';
    if(vertex->GetStoreFlag()) cqv='+';
    osv << cqv << vertex-> GetID() << '\0';
    std::string svid(osv.str());
    ostr << "       " << std::setw(6) << svid;
    //ostr << "      " << std::setw(4) << vertex-> GetID();
    ostr.unsetf(std::ios::fixed);
    ostr.setf(std::ios::scientific|std::ios::right|std::ios::showpoint);
    ostr << "- X(" << std::setw(9) << std::setprecision(2) 
	<< vertex-> GetPosition().x()/mm 
	<< "," << std::setw(9) << std::setprecision(2) 
	<< vertex-> GetPosition().y()/mm
	<< "," << std::setw(9) << std::setprecision(2) 
	<< vertex-> GetPosition().z()/mm 
	<< "," << std::setw(9) << std::setprecision(2) 
	<< vertex-> GetTime()/ns << ")";
    ostr.unsetf(std::ios::scientific);
    
    ostr << " @" << vertex-> GetVolumeName()
	<< "-" << vertex-> GetVolumeNumber();
  } 
  ostr << G4endl;
  
}

////////////////////////////////////////////////////////////////////
void G4MCTSimParticle::Print(std::ostream& ostr, G4bool qrevorder) const
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
