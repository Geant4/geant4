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
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4
//
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4HadronProcessStore.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int id = aRun->GetRunID();
  G4cout << "### Run " << id << " start" << G4endl;
  (HistoManager::GetPointer())->BeginOfRun();

#ifdef G4VIS_USE
  G4UImanager* UI = G4UImanager::GetUIpointer();

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run*)
{
  HistoManager* manager = HistoManager::GetPointer(); 
  G4cout << "RunAction: End of run actions are started" << G4endl;
  manager->EndOfRun();

#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
#endif

  if(manager->CrossSectionFlag()) {
    const G4Element* elm = manager->TargetElement();

    const G4ParticleDefinition* part[10];
    part[0] = G4Proton::Proton();
    part[1] = G4Neutron::Neutron();
    part[2] = G4PionPlus::PionPlus();
    part[3] = G4PionMinus::PionMinus();
    part[4] = G4KaonPlus::KaonPlus();
    part[5] = G4KaonMinus::KaonMinus();
    part[6] = G4AntiProton::AntiProton();
    part[7] = G4AntiNeutron::AntiNeutron();
    G4double emin = MeV;
    G4double emax = TeV;
    G4int nbins   = 60;
    G4double fac  = std::pow(emax/emin, 1.0/G4double(nbins));
    G4HadronProcessStore* store = G4HadronProcessStore::Instance();

    G4cout << "---------------------------------------------------------" 
	   << G4endl;
    for(G4int i=0; i<8; i++) {
      G4double e    = emin;
      G4cout << "---------------------------------------------------------" 
	     << G4endl;
      G4cout << "  N       e(GeV)      XsecInel(mb)     XsecEl(mb)   " 
	     << part[i]->GetParticleName() 
	     << G4endl;   
      G4cout << "---------------------------------------------------------" 
	     << G4endl;
      for(G4int j=0; j<nbins; j++) {
        G4double inel = store->GetInelasticCrossSectionPerAtom(part[i],e,elm)/millibarn; 
        G4double el = store->GetElasticCrossSectionPerAtom(part[i],e,elm)/millibarn; 
	G4cout << j << "  " << std::setw(8) << e/GeV << "  " 
	       << std::setw(8) << inel << "   " << std::setw(12) << el
	       << G4endl;
        manager->Fill(i*2 + 20, std::log10(e), inel);
        manager->Fill(i*2 + 21, std::log10(e), el);
	e *= fac;
      }
    }
    G4cout << "---------------------------------------------------------" 
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
