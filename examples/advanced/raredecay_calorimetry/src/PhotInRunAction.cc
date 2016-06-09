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
// $Id: PhotInRunAction.cc,v 1.5 2006/06/29 16:25:25 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

//#define debug

#include "PhotInRunAction.hh"

PhotInRunAction::PhotInRunAction() {}

PhotInRunAction::~PhotInRunAction() {}

G4Run* PhotInRunAction::GenerateRun() 
{
#ifdef debug
  G4cout<<"PhotInRunAction::GenerateRun: New Run is created"<<G4endl;
#endif
  return new PhotInRun; // @@ Who deletes it? (M.K.) What about previous Runs?
}

void PhotInRunAction::EndOfRunAction(const G4Run* aRun)
{
#ifdef debug
  G4cout<<"PhotInRunAction::EndOfRunAction is called"<<G4endl;
#endif
  const PhotInRun* theRun = (const PhotInRun*) aRun; // change the type
  // @@ Writing without any conditions (? M.K.)
  G4cout<<"##############################################################"<< G4endl;
  G4cout<<" Run Summary - Number of events : "<<theRun->GetNumberOfEvent()<<G4endl;
  G4cout<<"##############################################################"<< G4endl;
  G4double nEvt = (G4double)(theRun->GetNumberOfEvent());
  for(G4int i=0; i<PhotInNumSections; i++)
  {
    G4int ih1 = i+i;
    G4int ih2 = ih1+1;
    G4Region* region = G4RegionStore::GetInstance()->GetRegion(PhotInRegName[i]);
    G4ProductionCuts* cuts = region->GetProductionCuts();
    G4cout<<"Region "<<region->GetName()<<G4endl;
    G4cout<<"===================="<<G4endl;
    G4cout<<" Production thresholds :"<<G4endl;
    G4cout<<" -----------------------"<<G4endl;
    G4cout<<" gamma " <<G4BestUnit(cuts->GetProductionCut("gamma"),"Length")
          <<"    e- " <<G4BestUnit(cuts->GetProductionCut("e-"),"Length")
          <<"    e+ " <<G4BestUnit(cuts->GetProductionCut("e+"),"Length")<<G4endl;
    G4cout<<" Energy deposition in Absorber "
          <<G4BestUnit((theRun->GetTotalE(ih1))/nEvt,"Energy")<<", in Gap "
          <<G4BestUnit((theRun->GetTotalE(ih2))/nEvt,"Energy")<<G4endl;
    G4cout<<" Number of secondaries in the event :"<<G4endl;
    G4cout<<" ------------------------------------"<<G4endl;
    G4cout<<" gamma in Absorber " << (theRun->GetNGamma(ih1))/nEvt
          <<", in Gap      " << (theRun->GetNGamma(ih2))/nEvt<<G4endl;
    G4cout<<" e-    in Absorber " << (theRun->GetNElectron(ih1))/nEvt
          <<", in Gap " << (theRun->GetNElectron(ih2))/nEvt<<G4endl;
    G4cout<<" e+    in Absorber " << (theRun->GetNPositron(ih1))/nEvt
          <<", in Gap " << (theRun->GetNPositron(ih2))/nEvt<<G4endl;
    G4cout<<" Minimum kinetic energy of generated secondaries :"<<G4endl;
    G4cout<<" -------------------------------------------------"<<G4endl;
    G4cout<<" gamma in Absorber = "<<G4BestUnit(theRun->GetEMinGamma(ih1),"Energy")
          <<", in Gap = "<<G4BestUnit(theRun->GetEMinGamma(ih2),"Energy")<<G4endl;
    G4cout<<" e-    in Absorber = "<<G4BestUnit(theRun->GetEMinElectron(ih1),"Energy")
          <<", in Gap = " << G4BestUnit(theRun->GetEMinElectron(ih2),"Energy")<<G4endl;
    G4cout<<" e+    in Absorber = " << G4BestUnit(theRun->GetEMinPositron(ih1),"Energy")
          <<"    in Gap = " << G4BestUnit(theRun->GetEMinPositron(ih2),"Energy")<<G4endl;
    G4cout<<" Mean total track length of neutrons in the event : in Absorber = "
          <<G4BestUnit((theRun->GetTotalL(ih1))/nEvt,"Length")<<", in Gap = "
          <<G4BestUnit((theRun->GetTotalL(ih2))/nEvt,"Length")<<G4endl;
    G4cout<<" Mean number of steps of charged particles in the event : in Absorber ="
          <<theRun->GetNStep(ih1)/nEvt<<", in Gap = "<<theRun->GetNStep(ih2)/nEvt<<G4endl;
    G4cout<<"##############################################################"<<G4endl;
  }
}

