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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "G4Event.hh"
#include "Randomize.hh"

#include "MicrobeamEventAction.hh"
#include "MicrobeamRunAction.hh"
#include "MicrobeamHistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamEventAction::MicrobeamEventAction(MicrobeamRunAction* run,
MicrobeamHistoManager * his)
:Run(run),Histo(his),drawFlag("all"),printModulo(10000)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamEventAction::~MicrobeamEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamEventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  Run->SetNumEvent(evtNb);
  Run->SetDoseN(0);
  Run->SetDoseC(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamEventAction::EndOfEventAction(const G4Event* )
{  
 
// SAVE TOTAL ABSORBED DOSE IN PHANTOM

if (Run->GetDoseN()>0 || Run->GetDoseC()>0) 
{
	Histo->FillNtuple(3,0,Run->GetDoseN());
	Histo->FillNtuple(3,1,Run->GetDoseC());
        Histo->AddRowNtuple(3);			

        G4cout << "   ===> The incident alpha particle has reached the targeted cell :" << G4endl;
	G4cout << "   -----> total absorbed dose within Nucleus   is (Gy) = " << Run->GetDoseN() << G4endl;
	G4cout << "   -----> total absorbed dose within Cytoplasm is (Gy) = " << Run->GetDoseC() << G4endl;
	G4cout << G4endl;
}
else
{
	G4cout << "   ===> Sorry, the incident alpha particle has missed the targeted cell !" << G4endl;
	G4cout << G4endl;
}
}
