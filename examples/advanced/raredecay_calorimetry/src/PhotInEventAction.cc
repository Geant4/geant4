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
//
// $Id: PhotInEventAction.cc,v 1.6 2006/06/29 16:25:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

//#define debug

#include "PhotInEventAction.hh"

G4int PhotInEventAction::verboseLevel=0;

PhotInEventAction::PhotInEventAction()
{
#ifdef debug
  G4cout<<"PhotInEventAction::Constructor is called"<<G4endl;
#endif
  for(G4int i=0; i<PhotInDiNSections; i++) calorimeterCollID[i] = -1; // Reset CalorCollect
}

PhotInEventAction::~PhotInEventAction(){}

void PhotInEventAction::BeginOfEventAction(const G4Event*)
{
  for(G4int i=0; i<PhotInDiNSections; i++)
  {
    if(calorimeterCollID[i]==-1) // Name the unnamed CalorimeterCollections in the Begining
    {
      G4String colName=PhotInColNms[i];
#ifdef debug
      G4cout<<"PhotInEventAction::BeginOfEventAction:Col#"<<i<<",create="<<colName<<G4endl;
#endif
      calorimeterCollID[i] = G4SDManager::GetSDMpointer()->GetCollectionID(colName);
    }
  }
#ifdef debug
  G4cout<<"PhotInEventAction::BeginOfEventAction: End of Begin"<<G4endl;
#endif
}

void PhotInEventAction::EndOfEventAction(const G4Event* evt)
{
#ifdef debug
  G4cerr<<"PhotInEventAction::EndOfEventAction: is called"<<G4endl;
#endif
  if(verboseLevel==0) return; // Report of the Collection made in SD during this Event
  if(evt->GetEventID()>4 && (evt->GetEventID())%10>(verboseLevel-1)) return; // Only first

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return; // There is no collected information for this Event
#ifdef debug
  G4cerr<<"PhotInEventAction::EndOfEventAction:(HCE#0), Event="<<evt->GetEventID()<<G4endl;
#endif
  PhotInCalorHitsCollection* CHC = 0;
  for(G4int i=0; i<PhotInDiNSections; i++) // Make final sums for each collection and print
  {
    G4double totE=0.;
    G4double totL=0.;
    G4int    nStep=0;
    
    if (HCE) CHC = (PhotInCalorHitsCollection*)(HCE->GetHC(calorimeterCollID[i]));
    if (CHC)
    {
      G4int nHit = CHC->entries();
      for (G4int ii=0; ii<nHit; ii++)
      {
        totE  += (*CHC)[ii]->GetEDepos(); 
        totL  += (*CHC)[ii]->GetTrackL();
        nStep += (*CHC)[ii]->GetNSteps();
      }
    }
				G4cout<<PhotInColNms[i]<<" : "<<G4endl;
    G4cout<<" total energy deposition : "<<std::setw(7)<<G4BestUnit(totE,"Energy")<<G4endl;
    // @@ Change after modification of PhotInStackingAction
    //G4cout<<" number of particles generated :"<<G4endl
    //      <<"  gamma "<<PhotInStackingAction::GetNGamma(i) 
    //      <<"  e-    "<<PhotInStackingAction::GetNElectron(i) 
    //      <<"  e+    "<<PhotInStackingAction::GetNPositron(i)<< G4endl;
    // @@ Change after modification of PhotInStackingAction
    //G4cout<<" minimum kinetic energy of generated secondaries :"<<G4endl
    //      <<"  gamma "<<G4BestUnit(PhotInStackingAction::GetEMinGamma(i),"Energy") 
    //      <<"  e-    "<<G4BestUnit(PhotInStackingAction::GetEMinElectron(i),"Energy") 
    //      <<"  e+    "<<G4BestUnit(PhotInStackingAction::GetEMinPositron(i),"Energy")
    //      <<G4endl;
    G4cout<<" total track length of neutrons ="<<G4BestUnit(totL,"Length")<<" consists of "
          <<nStep<<" steps, meanStepLength="<<G4BestUnit(totL/nStep,"Length")<<G4endl;
    //@@ This can be written to the file !!
  }
#ifdef debug
  G4cerr<<"PhotInEventAction::EndOfEventAction: End of End"<<G4endl;
#endif
}  
