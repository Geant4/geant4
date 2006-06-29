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
// ====================================================================
//
//   H02EventAction.cc
//   $Id: H02EventAction.cc,v 1.4 2006-06-29 17:10:23 gunter Exp $
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
  const G4Event* ev = anEvent; ev=0;
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
void H02EventAction::EndOfEventAction(const G4Event*)
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
