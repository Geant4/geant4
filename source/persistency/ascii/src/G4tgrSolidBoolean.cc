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
// G4tgrSolidBoolean implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgrSolidBoolean.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrFileReader.hh"

// --------------------------------------------------------------------
G4tgrSolidBoolean::G4tgrSolidBoolean(const std::vector<G4String>& wl)
{
  // :SOLID/:VOLU VOLU UNION/SUBS/INTERS VOLU1 VOLU2 ROTM POSX POSY POSZ

  if(wl.size() != 9)
  {
    G4tgrUtils::DumpVS(wl, "G4tgrSolidBoolean::G4tgrSolidBoolean()");
    G4Exception("G4tgrSolidBoolean::G4tgrSolidBoolean()", "InvalidInput",
                FatalException, "Line read with less or more than 9 words.");
  }

  //---------- Set name
  theName = G4tgrUtils::GetString(wl[1]);

  G4tgrVolumeMgr* volmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrSolid* sol1 = volmgr->FindSolid(G4tgrUtils::GetString(wl[3]));
  if(sol1 == nullptr)
  {
    sol1 = volmgr->FindVolume(G4tgrUtils::GetString(wl[3]), 1)->GetSolid();
  }
  const G4tgrSolid* sol2 = volmgr->FindSolid(G4tgrUtils::GetString(wl[4]));
  if(sol2 == nullptr)
  {
    sol2 = volmgr->FindVolume(G4tgrUtils::GetString(wl[4]), 1)->GetSolid();
  }
  theSolids.push_back(sol1);
  theSolids.push_back(sol2);

  //---------- Set relative placement and rotation matrix
  theRelativeRotMatName = G4tgrUtils::GetString(wl[5]);
  theRelativePlace =
    G4ThreeVector(G4tgrUtils::GetDouble(wl[6]), G4tgrUtils::GetDouble(wl[7]),
                  G4tgrUtils::GetDouble(wl[8]));
  //---------- Set solid type
  G4String wl2 = wl[2];
  for(G4int ii = 0; ii < (G4int)wl2.length(); ++ii)
  {
    wl2[ii] = (char)std::toupper(wl2[ii]);
  }
  theType = "Boolean_" + wl2;

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif

  G4tgrVolumeMgr::GetInstance()->RegisterMe(this);
}

// --------------------------------------------------------------------
G4tgrSolidBoolean::~G4tgrSolidBoolean()
{
}

// --------------------------------------------------------------------
const G4String& G4tgrSolidBoolean::GetRelativeRotMatName() const
{
  return theRelativeRotMatName;
}

// --------------------------------------------------------------------
G4ThreeVector G4tgrSolidBoolean::GetRelativePlace() const
{
  return theRelativePlace;
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrSolidBoolean& sol)
{
  os << "G4tgrSolidBoolean= " << sol.theName << " of type " << sol.theType
     << " PARAMS: ";
  if(sol.theSolidParams.size() != 0)
  {
    std::vector<G4double> solpar = *(sol.theSolidParams[0]);
    for(std::size_t ii = 0; ii < solpar.size(); ++ii)
    {
      os << solpar[ii] << " ";
    }
  }
  os << G4endl;

  return os;
}
