// $Id: G4MCTSimVertex.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
// ====================================================================
//
//   G4MCTSimVertex.cc
//
// ====================================================================
#include <strstream>
#include <iomanip>
#include "CLHEP/Units/SystemOfUnits.h"
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
	       std::string vname, int ncopy, std::string pname)
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
void G4MCTSimVertex::Print(std::ostream& ostr) const
//////////////////////////////////////////////////
{
  char st[20];
  std::ostrstream os(st,20);
  char cq=' ';
  if(storeFlag) cq='+';
  os << cq << id << '\0';
  std::string sid(st);
  
  ostr.unsetf(std::ios::fixed);
  ostr.setf(std::ios::scientific|std::ios::right|std::ios::showpoint);
  //ostr << std::setw(4) << id;
  ostr << std::setw(6) << sid;
  ostr << " : X(" << std::setw(9) << std::setprecision(2) 
       << position.x()/mm 
       << "," << std::setw(9) << std::setprecision(2) 
       << position.y()/mm
       << "," << std::setw(9) << std::setprecision(2) 
       << position.z()/mm 
       << "," << std::setw(9) << std::setprecision(2) 
       << time/ns << ")";
  ostr.unsetf(std::ios::scientific);
  ostr << "@" << volumeName
       << "-" << volumeNumber
       << "%" << creatorProcessName
       << std::endl;
  
  ostr << "      " << std::setw(4) << inParticleTrackID << "-> ";
  int np= outParticleTrackIDList.size();
  for (int i=0; i<np; i++) ostr << outParticleTrackIDList[i] << ", ";
  ostr << std::endl;
}

