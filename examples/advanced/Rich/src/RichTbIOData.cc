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
// Rich advanced example for Geant4
// RichTbIOData.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbIOData.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "RichTbHit.hh"
#include <fstream>

RichTbIOData::RichTbIOData(RichTbRunConfig* RConfig ) 
  :OutputDataFS((RConfig -> getOutputFileName()).c_str())
{

   aOutFileString=RConfig -> getOutputFileName();
  
   
}
RichTbIOData::~RichTbIOData() { }

void  RichTbIOData::WriteOutEventHeaderData( const G4Event* ) { }

void  RichTbIOData::WriteOutHitData( const G4Event* evt ) { 

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int RichTbCID = SDman->GetCollectionID(colNam="RichTbHitsCollection");
  if(RichTbCID<0) return;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  RichTbHitsCollection* RHC = NULL;
  if(HCE)
  {
    RHC = (RichTbHitsCollection*)(HCE->GetHC(RichTbCID));
  }
  if(RHC)
  {
    G4int n_hit = RHC->entries();
    G4cout << "     " << n_hit << "Hits being written out "<<G4endl;
    // OutputDataFS<<n_hit<<G4endl;
    for (G4int ih=0; ih<n_hit; ih++ ) {
      RichTbHit* aHit = (*RHC)[ih];
      G4int aHitHpdNum = aHit -> GetCurHpdNum();
      G4int aHitSectNum = aHit -> GetCurSectNum();
      G4int aHitPixelNum = aHit -> GetCurPixNum();
      G4ThreeVector aHitPosPC = aHit -> GetPosPC();

      OutputDataFS<<"  "<< aHitHpdNum<<"  "<< aHitSectNum
                  <<"  "<<aHitPixelNum
                  <<"  "<< aHitPosPC.x()<<"  "<< aHitPosPC.y()
                  <<"  "<< aHitPosPC.z()<<G4endl;


    }
  }


}

