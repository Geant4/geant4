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
// $Id: ExN07RunAction.cc,v 1.1 2003-03-10 01:43:37 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN07RunAction.hh"
#include "ExN07Run.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

ExN07RunAction::ExN07RunAction()
{;}

ExN07RunAction::~ExN07RunAction()
{;}

G4Run* ExN07RunAction::GenerateRun()
{ return new ExN07Run; }

void ExN07RunAction::EndOfRunAction(const G4Run* aRun)
{
  static G4String regName[3] = {"Calor-A","Calor-B","Calor-C"};

  const ExN07Run* theRun = (const ExN07Run*)aRun;
  
  G4cout
    << "############################################################" << G4endl;
  G4cout 
    << " Run Summary - Number of events : " << theRun->GetNumberOfEvent() 
    << G4endl;
  G4cout
    << "############################################################" << G4endl;

  G4double nEvt = (G4double)(theRun->GetNumberOfEvent());
  for(size_t i=0;i<3;i++)
  {
    size_t ih1 = 2*i;
    size_t ih2 = 2*i+1;

    G4Region* region = G4RegionStore::GetInstance()->GetRegion(regName[i]);
    G4ProductionCuts* cuts = region->GetProductionCuts();
    G4cout << "Region : " << region->GetName() << G4endl;
    G4cout << "  Production thresholds :" 
           << " gamma " << G4BestUnit(cuts->GetProductionCut("gamma"),"Length")
           << "    e- " << G4BestUnit(cuts->GetProductionCut("e-"),"Length")
           << "    e+ " << G4BestUnit(cuts->GetProductionCut("e+"),"Length")
           << G4endl;
    G4cout << "  Average energy deposition :"
           << " Absorber " << G4BestUnit((theRun->GetTotalE(ih1))/nEvt,"Energy")
           << "      Gap " << G4BestUnit((theRun->GetTotalE(ih2))/nEvt,"Energy")
           << G4endl;
    G4cout << "  Average track length of e+/e- :"
           << " Absorber " << G4BestUnit((theRun->GetTotalL(ih1))/nEvt,"Length")
           << "      Gap " << G4BestUnit((theRun->GetTotalL(ih2))/nEvt,"Length")
           << G4endl;
    G4cout << "  Average number of steps of e+/e- :"
           << " Absorber " << (theRun->GetNStep(ih1))/nEvt
           << "      Gap " << (theRun->GetNStep(ih2))/nEvt
           << G4endl;
    G4cout
    << "############################################################" << G4endl;
  }
}

