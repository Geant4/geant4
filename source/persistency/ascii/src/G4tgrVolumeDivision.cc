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
// G4tgrVolumeDivision implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgrVolumeDivision.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlace.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrMessenger.hh"

// --------------------------------------------------------------------
G4tgrVolumeDivision::~G4tgrVolumeDivision()
{
}

// --------------------------------------------------------------------
G4tgrVolumeDivision::G4tgrVolumeDivision(const std::vector<G4String>& wl)
{
  // wl: NAME PARENT  MATERIAL AXIS STEP/NDIV OFFSET

  G4tgrUtils::CheckWLsize(wl, 6, WLSIZE_GE,
                          "G4tgrVolumeDivision::G4tgrVolumeDivision");
  G4tgrUtils::CheckWLsize(wl, 8, WLSIZE_LE,
                          "G4tgrVolumeDivision::G4tgrVolumeDivision");

  theType = "VOLDivision";

  // :DIV  NAME PARENT MATERIAL AXIS STEP/NDIV OFFSET

  //---------- set name
  theName = G4tgrUtils::GetString(wl[1]);

  //---------- set the pointer to the parent DU
  G4String parentName = G4tgrUtils::GetString(wl[2]);
  G4tgrVolumeMgr::GetInstance()->FindVolume(parentName, 1);  // check existance

  //---------- initialize G4tgrPlace
  thePlaceDiv = new G4tgrPlaceDivRep();
  thePlaceDiv->SetParentName(parentName);
  thePlaceDiv->SetType("PlaceDivision");
  thePlaceDiv->SetVolume(this);

  //---------- set material name
  theMaterialName = G4tgrUtils::GetString(wl[3]);

  //----- set axis of replica
  thePlaceDiv->SetAxis(thePlaceDiv->BuildAxis(G4tgrUtils::GetString(wl[4])));

  //------ register parent - child
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild(parentName, thePlaceDiv);
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 3)
  {
    G4cout << " G4tgrVolumeDivision::G4tgrVolumeDivision() -"
           << " Replica register parent - child " << G4endl;
  }
#endif

  //---------- set if division is given by number of divisions of by width
  G4String wl0 = wl[0];
  for(G4int ii = 0; ii < (G4int)wl0.length(); ++ii)
  {
    wl0[ii] = (char)std::toupper(wl0[ii]);
  }

  if(wl0 == ":DIV_NDIV")
  {
    thePlaceDiv->SetDivType(DivByNdiv);
    thePlaceDiv->SetNDiv(G4tgrUtils::GetInt(wl[5]));
    if(wl.size() == 7)
    {
      thePlaceDiv->SetOffset(G4tgrUtils::GetDouble(wl[6]) * mm);
    }
  }
  else if(wl0 == ":DIV_WIDTH")
  {
    thePlaceDiv->SetDivType(DivByWidth);
    thePlaceDiv->SetWidth(G4tgrUtils::GetDouble(wl[5]) * mm);
    if(wl.size() == 7)
    {
      thePlaceDiv->SetOffset(G4tgrUtils::GetDouble(wl[6]) * mm);
    }
  }
  else if(wl0 == ":DIV_NDIV_WIDTH")
  {
    thePlaceDiv->SetDivType(DivByNdivAndWidth);
    thePlaceDiv->SetNDiv(G4tgrUtils::GetInt(wl[5]));
    thePlaceDiv->SetWidth(G4tgrUtils::GetDouble(wl[6]) * mm);
    if(wl.size() == 8)
    {
      thePlaceDiv->SetOffset(G4tgrUtils::GetDouble(wl[7]) * mm);
    }
  }
  else
  {
    G4String ErrMessage = "Division type not supported, sorry... " + wl[0];
    G4Exception("G4tgrVolumeDivision::G4tgrVolumeDivision()", "NotImplemented",
                FatalException, ErrMessage);
  }

  theVisibility = 1;
  theRGBColour  = new G4double[3];
  for(std::size_t ii = 0; ii < 3; ++ii)
  {
    theRGBColour[ii] = -1.;
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif

  theSolid = nullptr;
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrVolumeDivision& obj)
{
  os << "G4tgrVolumeDivision= " << obj.theName
     << " Placement= " << *(obj.thePlaceDiv) << G4endl;

  return os;
}
