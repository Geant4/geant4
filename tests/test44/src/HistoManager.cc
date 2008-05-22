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
// $Id: HistoManager.cc,v 1.4 2008-05-22 15:54:51 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
//
// Author:      V.Ivanchenko 30/01/01
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
// 16.11.2006 Add beamFlag (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "Histo.hh"
#include "G4Track.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager* HistoManager::GetPointer()
{
  if(!fManager) {
    static HistoManager manager;
    fManager = &manager;
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::HistoManager()
{
  verbose=  0;
  nSlices   = 300;
  nBinsE    = 100;
  nHisto    = 1;
  length    = 300.*mm;
  edepMax   = 1.0*GeV;
  beamFlag  = true;
  material  = 0;
  elm       = 0;
  histo     = new Histo();
  neutron   = G4Neutron::Neutron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::bookHisto()
{
  histo->add1D("1","Energy deposition (MeV/mm/event) in the target",
	       nSlices,0.0,length/mm,MeV/mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  absZ0       = -0.5*length;
  n_evt       = 0;
  n_step      = 0;

  edepSum     = 0.0;
  edepSum2    = 0.0;

  histo->setVerbose(verbose);
  bookHisto();
  histo->book();

  if(verbose > 0) {
    G4cout << "HistoManager: Histograms are booked and run has been started"
           <<G4endl<<"  BeginOfRun (After histo->book)"<< G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{

  G4cout << "HistoManager: End of run actions are started" << G4endl;

  // Average values
  G4cout<<"========================================================"<<G4endl;

  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;

  G4double xs = x*(G4double)n_step;

  edepSum  *= x;
  edepSum2 *= x;
  edepSum2 -= edepSum*edepSum;
  if(edepSum2 > 0.0) edepSum2 = std::sqrt(edepSum2);
  else               edepSum2 = 0.0;

  G4cout                         << "Beam particle                        "
				 << primaryDef->GetParticleName() <<G4endl;
  G4cout                         << "Beam Energy(MeV)                     " 
				 << primaryKineticEnergy/MeV <<G4endl;
  G4cout                         << "Number of events                     " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit (MeV)         " << edepSum/MeV 
	 << "   RMS(MeV) " << edepSum2/MeV << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps              " << xs << G4endl;
  G4cout<<"========================================================"<<G4endl;
  G4cout<<G4endl;

  // normalise histograms
  for(G4int i=0; i<nHisto; i++) {
    histo->scale(i,x);
  }

  histo->print(0);
//  histo->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfEvent()
{
  edepEvt = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfEvent()
{
  edepSum  += edepEvt;
  edepSum2 += edepEvt*edepEvt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::ScoreNewTrack(const G4Track* track)
{
  const G4ParticleDefinition* pd = track->GetDefinition();
  G4String name = pd->GetParticleName();
  G4double e = track->GetKineticEnergy();

  // Primary track
  if(0 == track->GetParentID()) {

    n_evt++;
    primaryKineticEnergy = e;
    primaryDef = pd;
    G4ThreeVector dir = track->GetMomentumDirection();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddTargetStep(const G4Step* step)
{
  n_step++;
  G4double edep = step->GetTotalEnergyDeposit();
  if(edep >= DBL_MIN) { 
    const G4Track* track = step->GetTrack();
    currentDef = track->GetDefinition(); 
    currentKinEnergy = track->GetKineticEnergy();

    G4ThreeVector pos = 
      (step->GetPreStepPoint()->GetPosition() +
       step->GetPostStepPoint()->GetPosition())*0.5;

    G4double z = pos.z() - absZ0;

    // scoring
    edepEvt += edep;
    histo->fill(0,z,edep);

    if(1 < verbose) 
      G4cout << "HistoManager::AddEnergy: e(keV)= " << edep/keV
	     << "; z(mm)= " << z/mm
	     << "; step(mm)= " << step->GetStepLength()/mm
	     << " by " << currentDef->GetParticleName()
	     << " E(MeV)= " << currentKinEnergy/MeV
	     << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  verbose = val; 
  histo->setVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetTargetMaterial(const G4Material* mat)         
{
  if(mat) {
    material = mat;
    elm = (*(material->GetElementVector()))[0];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::Fill(G4int id, G4double x, G4double w)
{
  histo->fill(id, x, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

