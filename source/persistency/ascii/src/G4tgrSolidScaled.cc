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
// G4tgrSolidScaled implementation
//
// Author: P.Heidary, AEOI - November 2021
// --------------------------------------------------------------------

#include "G4tgrSolidScaled.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrMessenger.hh"

// --------------------------------------------------------------------
G4tgrSolidScaled::G4tgrSolidScaled(const std::vector<G4String>& wl)
{
  // :SOLID/:VOLU SOL1 SCALED SOL0 SCLAEX SCALEY SCALEZ
  if(wl.size() != 7)
  {
    G4tgrUtils::DumpVS(wl, "G4tgrSolidScaled::G4tgrSolidScaled()");
    G4Exception("G4tgrSolidScaled::G4tgrSolidScaled()", "InvalidInput",
                FatalException, "Line read with less or more than 7 words.");
  }

  //---------- Set name
  theName = G4tgrUtils::GetString(wl[1]);

  G4tgrVolumeMgr* volmgr = G4tgrVolumeMgr::GetInstance();
  origSolid = volmgr->FindSolid(G4tgrUtils::GetString(wl[3]));
  if(origSolid == nullptr)
  {
    origSolid = volmgr->FindVolume(G4tgrUtils::GetString(wl[3]), 1)->GetSolid();
  }

  //---------- Set scale vector
  scale3d = G4Scale3D(G4tgrUtils::GetDouble(wl[4]), G4tgrUtils::GetDouble(wl[5]), 
    G4tgrUtils::GetDouble(wl[6]));

  //---------- Set solid type
  G4String wl2 = wl[2];
  for(G4int ii = 0; ii < (G4int)wl2.length(); ++ii)
  {
    wl2[ii] = (char)std::toupper(wl2[ii]);
  }
  theType = wl2;

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif

  G4tgrVolumeMgr::GetInstance()->RegisterMe(this);
}

// --------------------------------------------------------------------
G4tgrSolidScaled::~G4tgrSolidScaled()
{
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrSolidScaled& sol)
{
  os << "G4tgrSolidScaled= " << sol.theName << " of type " << sol.theType
     << " original solid: " << sol.origSolid->GetName() << " Scale x: " << 
     sol.scale3d.xx() << " Scale y: "     << sol.scale3d.yy() << 
     " Scale z: " << sol.scale3d.zz() << G4endl;

  return os;
}
