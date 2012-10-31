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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"

#include "IORTEventAction.hh"
#include "IORTDetectorHit.hh"
#include "IORTDetectorSD.hh"
#include "IORTDetectorConstruction.hh"
#include "IORTMatrix.hh"
#include "IORTEventActionMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
IORTEventAction::IORTEventAction() :
  drawFlag("all" ),printModulo(1000), pointerEventMessenger(0)
{ 
  hitsCollectionID = -1;
  pointerEventMessenger = new IORTEventActionMessenger(this);
}

/////////////////////////////////////////////////////////////////////////////
IORTEventAction::~IORTEventAction()
{
 delete pointerEventMessenger;
}

/////////////////////////////////////////////////////////////////////////////
void IORTEventAction::BeginOfEventAction(const G4Event* evt)
{ 
  G4int evtNb = evt->GetEventID();
  
  //printing survey
  if (evtNb%printModulo == 0)
     G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
   
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
  if(hitsCollectionID == -1)
    hitsCollectionID = pSDManager -> GetCollectionID("IORTDetectorHitsCollection");
  
}

/////////////////////////////////////////////////////////////////////////////
void IORTEventAction::EndOfEventAction(const G4Event* evt)
{ 
  if(hitsCollectionID < 0)
  return;
  G4HCofThisEvent* HCE = evt -> GetHCofThisEvent();

  // Clear voxels hit list 
  IORTMatrix* matrix = IORTMatrix::GetInstance();
  if (matrix) matrix -> ClearHitTrack(); 

  if(HCE)
  {
    IORTDetectorHitsCollection* CHC = (IORTDetectorHitsCollection*)(HCE -> GetHC(hitsCollectionID));
    if(CHC)
     {
       if(matrix)
	  { 
	      // Fill the matrix with the information: voxel and associated energy deposit 
          // in the detector at the end of the event

	  G4int HitCount = CHC -> entries();
	  for (G4int h=0; h<HitCount; h++)
	    {
	      G4int i = ((*CHC)[h]) -> GetXID();
	      G4int j = ((*CHC)[h]) -> GetYID();
	      G4int k = ((*CHC)[h]) -> GetZID();
              G4double energyDeposit = ((*CHC)[h]) -> GetEdep();
              matrix -> Fill(i, j, k, energyDeposit/MeV);              
	    }
	  }
    }
  }
}

