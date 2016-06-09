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
//   G4MCTSimVertex.cc
//
// ====================================================================

#include <sstream>
#include <iomanip>

#include "G4MCTSimVertex.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4MCTSimParticle.hh"

// ====================================================================
//
// class description
//
// ====================================================================

/////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex()  
  : inParticleTrackID(0),id(-1), time(0.),
    volumeName(""), volumeNumber(-1),
    creatorProcessName("none"),
    storeFlag(false)
/////////////////////////////////////
{
}

/////////////////////////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex(const G4ThreeVector& x, double t)
  : inParticleTrackID(0), id(-1), position(x), time(t),
    volumeName(""), volumeNumber(-1),
    creatorProcessName("none"), storeFlag(false)
/////////////////////////////////////////////////////////
{
}

///////////////////////////////////////////////////////////////
G4MCTSimVertex::G4MCTSimVertex(const G4ThreeVector& x, double t,
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
  std::ostringstream os;
  char cq=' ';
  if(storeFlag) cq='+';
  os << cq << id << '\0';
  std::string sid(os.str());
  
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
       << G4endl;
  
  ostr << "      " << std::setw(4) << inParticleTrackID << "-> ";
  size_t np= outParticleTrackIDList.size();
  for (size_t i=0; i<np; i++) ostr << outParticleTrackIDList[i] << ", ";
  ostr << G4endl;
}

