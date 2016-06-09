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
//   G4MCTSimVertex.cc
//
// ====================================================================
#include "g4std/strstream"
#include "g4std/iomanip"
#include "G4ios.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4MCTSimParticle.hh"
#include "G4MCTSimVertex.hh"

// ====================================================================
//
// class description
//
// ====================================================================

/////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex()  
  : inParticleTrackID(0),id(-1), 
    volumeName(""), volumeNumber(-1),
    creatorProcessName("none"),
    storeFlag(false)
/////////////////////////////////////
{
}

/////////////////////////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex(const Hep3Vector& x, double t)
  : inParticleTrackID(0), id(-1), position(x), time(t),
    volumeName(""), volumeNumber(-1),
    creatorProcessName("none"), storeFlag(false)
/////////////////////////////////////////////////////////
{
}

///////////////////////////////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex(const Hep3Vector& x, double t,
	       G4std::string vname, int ncopy, G4std::string pname)
  : inParticleTrackID(0), id(-1), position(x), time(t),
    volumeName(vname), volumeNumber(ncopy),
    creatorProcessName(pname), storeFlag(false)
///////////////////////////////////////////////////////////////
{
}

/////////////////////////////
G4MCTSimVertex::~G4MCTSimVertex()
/////////////////////////////
{
  outParticleTrackIDList.clear();
}

//////////////////////////////////////////////////
void G4MCTSimVertex::Print(G4std::ostream& ostr) const
//////////////////////////////////////////////////
{
  char st[20];
  G4std::ostrstream os(st,20);
  char cq=' ';
  if(storeFlag) cq='+';
  os << cq << id << '\0';
  G4std::string sid(st);
  
  ostr.unsetf(G4std::ios::fixed);
  ostr.setf(G4std::ios::scientific|G4std::ios::right|G4std::ios::showpoint);
  //ostr << G4std::setw(4) << id;
  ostr << G4std::setw(6) << sid;
  ostr << " : X(" << G4std::setw(9) << G4std::setprecision(2) 
       << position.x()/mm 
       << "," << G4std::setw(9) << G4std::setprecision(2) 
       << position.y()/mm
       << "," << G4std::setw(9) << G4std::setprecision(2) 
       << position.z()/mm 
       << "," << G4std::setw(9) << G4std::setprecision(2) 
       << time/ns << ")";
  ostr.unsetf(G4std::ios::scientific);
  ostr << "@" << volumeName
       << "-" << volumeNumber
       << "%" << creatorProcessName
       << G4endl;
  
  ostr << "      " << G4std::setw(4) << inParticleTrackID << "-> ";
  size_t np= outParticleTrackIDList.size();
  for (size_t i=0; i<np; i++) ostr << outParticleTrackIDList[i] << ", ";
  ostr << G4endl;
}

