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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"

#include "HadrontherapyEventAction.hh"
#include "HadrontherapyDetectorHit.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyEventActionMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapyEventAction::HadrontherapyEventAction() :
drawFlag("all" ),printModulo(10), pointerEventMessenger(0)
{
    hitsCollectionID = -1;
    pointerEventMessenger = new HadrontherapyEventActionMessenger(this);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyEventAction::~HadrontherapyEventAction()
{
    delete pointerEventMessenger;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyEventAction::BeginOfEventAction(const G4Event*)
{
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(hitsCollectionID == -1)
        hitsCollectionID = pSDManager -> GetCollectionID("HadrontherapyDetectorHitsCollection");
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyEventAction::EndOfEventAction(const G4Event* evt)
{
    if(hitsCollectionID < 0)
        return;
    G4HCofThisEvent* HCE = evt -> GetHCofThisEvent();
    
    // Clear voxels hit list
    HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
    if (matrix) matrix -> ClearHitTrack();
    
    if(HCE)
    {
        HadrontherapyDetectorHitsCollection* CHC = (HadrontherapyDetectorHitsCollection*)(HCE -> GetHC(hitsCollectionID));
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
