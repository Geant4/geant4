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

void  RichTbIOData::WriteOutEventHeaderData( const G4Event* evt ) { }

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

