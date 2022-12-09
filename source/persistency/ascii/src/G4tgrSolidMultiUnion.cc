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
// G4tgrSolidMultiUnion implementation
//
// Author: P.Heidary, AEOI - November 2021
// --------------------------------------------------------------------

#include "G4tgrSolidMultiUnion.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbRotationMatrixMgr.hh"

// --------------------------------------------------------------------
G4tgrSolidMultiUnion::G4tgrSolidMultiUnion(const std::vector<G4String>& wl)
{
  // :SOLID/:VOLU SOL MULTIUNION NSOL SOL1 ROTM1 POSX1 POSY1 POSZ1
  //                              SOL2 ROTM2 POSX2 POSY2 POSZ2 ...

  //---------- Set name
  theName = G4tgrUtils::GetString(wl[1]);
  nSolid  = G4tgrUtils::GetInt(wl[3]);
  G4tgrVolumeMgr* volmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrSolid* sol1;

  if(G4int(wl.size()) != 4 + 5*nSolid)
  {
    G4String Err1 = "Solid type MULTIUNION with ";
    G4String Err2 = std::to_string(nSolid);
    G4String Err3 = " Solids, should have ";
    G4String Err4 = std::to_string(4 + 5*nSolid);
    G4String Err5 = " parameters.";
    G4String Err6 = " It has "+G4UIcommand::ConvertToString(int(wl.size()));
    
    G4String ErrMessage = Err1 + Err2 + Err3 + Err4 + Err5 + Err6;

    G4tgrUtils::DumpVS(wl, "G4tgrSolidMultiUnion::G4tgrSolidMultiUnion()");
    G4Exception("G4tgrSolidMultiUnion::G4tgrSolidMultiUnion()", "InvalidInput",
                FatalException, ErrMessage);
  }

  for (G4int i = 4; i< nSolid*5; i+=5)
  {
    sol1 = volmgr->FindSolid(G4tgrUtils::GetString(wl[i]));
    if(sol1 == nullptr)
    {
      sol1 = volmgr->FindVolume(G4tgrUtils::GetString(wl[i]), 1)->GetSolid();
    }
    theSolids.push_back(sol1);
    
    //---------- Populate transformation matrix
    theRotMatName = G4tgrUtils::GetString(wl[i + 1]);
  theRotMat =
    G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(theRotMatName);
  thePosition =
    G4ThreeVector(G4tgrUtils::GetDouble(wl[i+2]), 
                  G4tgrUtils::GetDouble(wl[i+3]),
                  G4tgrUtils::GetDouble(wl[i+4]));
  tr1 = G4Transform3D(*theRotMat, thePosition);
  theTransformations.push_back(tr1);
  }

  // ---------- Set solid type
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
G4tgrSolidMultiUnion::~G4tgrSolidMultiUnion()
{
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrSolidMultiUnion& sol)
{
  os << "G4tgrSolidMultiUnion= " << sol.theName << " of type " << sol.theType
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
