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

// ====================================================================
//
//   H02EventAction.cc
//   $Id: H02EventAction.cc,v 1.1 2002-05-28 14:15:46 murakami Exp $
//
// ====================================================================
#include "H02EventAction.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "H02MuonSD.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////
H02EventAction::H02EventAction()
//////////////////////////////
{
}

///////////////////////////////
H02EventAction::~H02EventAction()
///////////////////////////////
{
}


//////////////////////////////////////////////////////////////
void H02EventAction::BeginOfEventAction(const G4Event* anEvent)
//////////////////////////////////////////////////////////////
{
#ifdef DEBUG_HEPMC
  // printout primary information
  G4cout << "Print out primary information" << G4endl;
  G4int nVtx= anEvent-> GetNumberOfPrimaryVertex();
  G4int i;
  for(i=0; i< nVtx; i++) {
    const G4PrimaryVertex* primaryVertex= anEvent-> GetPrimaryVertex(i);
    primaryVertex-> Print();  
  }
#endif
}

////////////////////////////////////////////////////////////
void H02EventAction::EndOfEventAction(const G4Event* anEvent)
////////////////////////////////////////////////////////////
{
  G4cout << " Print out hit information" << G4endl;
  G4SDManager* SDManager= G4SDManager::GetSDMpointer();
  H02MuonSD* muonSD= 
    (H02MuonSD*)SDManager-> FindSensitiveDetector("/mydet/muon");
  muonSD-> PrintAll();
  muonSD-> DrawAll();
  G4cout << G4endl;
}
