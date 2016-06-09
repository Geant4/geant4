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
//
// $Id: ExN07EventAction.cc,v 1.5 2003/06/16 16:50:07 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//

#include "ExN07EventAction.hh"

#include "ExN07CalorHit.hh"
#include "ExN07StackingAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

G4int ExN07EventAction::verboseLevel=0;

ExN07EventAction::ExN07EventAction()
{
  for(size_t i=0;i<6;i++)
  { calorimeterCollID[i] = -1; }
}

ExN07EventAction::~ExN07EventAction()
{;}

void ExN07EventAction::BeginOfEventAction(const G4Event*)
{
  for(size_t i=0;i<6;i++)
  {
    if(calorimeterCollID[i]==-1)
    {
      G4String colName;
      switch(i)
      {
        case 0:
          colName = "CalorSD-A/AbsCollection"; break;
        case 1:
          colName = "CalorSD-A/GapCollection"; break;
        case 2:
          colName = "CalorSD-B/AbsCollection"; break;
        case 3:
          colName = "CalorSD-B/GapCollection"; break;
        case 4:
          colName = "CalorSD-C/AbsCollection"; break;
        case 5:
          colName = "CalorSD-C/GapCollection"; break;
      }
      calorimeterCollID[i] = 
        G4SDManager::GetSDMpointer()->GetCollectionID(colName);
    }
  }
}

void ExN07EventAction::EndOfEventAction(const G4Event* evt)
{
  if(verboseLevel==0) return;
  if(evt->GetEventID()>4 && (evt->GetEventID())%10>(verboseLevel-1)) return;

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  G4cout << G4endl << "Event : " << evt->GetEventID() << G4endl;

  ExN07CalorHitsCollection* CHC = 0;
  for(size_t i=0;i<6;i++)
  {
    G4double totE=0.;
    G4double totL=0.;
    G4int nStep=0;
    
    if (HCE) CHC = (ExN07CalorHitsCollection*)(HCE->GetHC(calorimeterCollID[i]));
    if (CHC)
    {
      G4int nHit = CHC->entries();
      for (G4int ii=0;ii<nHit;ii++)
      {
        totE += (*CHC)[ii]->GetEdep(); 
        totL += (*CHC)[ii]->GetTrak();
        nStep += (*CHC)[ii]->GetNStep();
      }
    }
   
    switch(i)
    {
      case 0:
        G4cout << "Calor-A : Absorber" << G4endl; break;
      case 1:
        G4cout << "Calor-A : SensitiveGap" << G4endl; break;
      case 2:
        G4cout << "Calor-B : Absorber" << G4endl; break;
      case 3:
        G4cout << "Calor-B : SensitiveGap" << G4endl; break;
      case 4:
        G4cout << "Calor-C : Absorber" << G4endl; break;
      case 5:
        G4cout << "Calor-C : SensitiveGap" << G4endl; break;
    }
    G4cout
       << "  total energy deposition : " << std::setw(7)
       << G4BestUnit(totE,"Energy") << G4endl;
    G4cout
       << "  number of particles generated :" << G4endl
       << "    gamma " << ExN07StackingAction::GetNGamma(i) 
       << "    e- " << ExN07StackingAction::GetNElectron(i) 
       << "    e+ " << ExN07StackingAction::GetNPositron(i) << G4endl;
    G4cout
       << "  minimum kinetic energy of generated secondaries :" << G4endl << std::setw(7)
       << "    gamma " << G4BestUnit(ExN07StackingAction::GetEMinGamma(i),"Energy") 
       << "    e- " << G4BestUnit(ExN07StackingAction::GetEMinElectron(i),"Energy") 
       << "    e+ " << G4BestUnit(ExN07StackingAction::GetEMinPositron(i),"Energy")
       << G4endl;
    G4cout
       << "  total track length of e+/e- : " << std::setw(7)
       << G4BestUnit(totL,"Length") << G4endl;
    G4cout
       << "  number of steps of e+/e- : " << nStep
       << G4endl;
  }
}  

