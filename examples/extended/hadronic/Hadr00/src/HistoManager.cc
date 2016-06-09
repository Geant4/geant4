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
/// \file hadronic/Hadr00/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id$
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
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4Neutron.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4NistManager.hh"
#include "G4HadronicProcessStore.hh"

#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include "G4StableIsotopes.hh"
#include "G4SystemOfUnits.hh"

#include "Histo.hh"

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
  fHisto     = new Histo();
  fNeutron   = G4Neutron::Neutron();
  fVerbose   = 1;

  fParticleName  = "proton";
  fElementName   = "Al";

  fMinKinEnergy  = 0.1*MeV;
  fMaxKinEnergy  = 10*TeV;
  fMinMomentum   = 1*MeV;
  fMaxMomentum   = 10*TeV;

  fBinsE    = 800;
  fBinsP    = 700;

  fHisto->SetVerbose(fVerbose);

  fIsInitialised = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete fHisto;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  G4double p1 = std::log10(fMinMomentum/GeV);
  G4double p2 = std::log10(fMaxMomentum/GeV);
  G4double e1 = std::log10(fMinKinEnergy/MeV);
  G4double e2 = std::log10(fMaxKinEnergy/MeV);

  //G4cout<<"e1= "<<e1<<" e2= "<<e2<<" p1= "<<p1<<" p2= "<<p2<<G4endl;

  if(fIsInitialised) {
    fHisto->SetHisto1D(0,fBinsP,p1,p2,1.0);
    fHisto->SetHisto1D(1,fBinsE,e1,e2,1.0);
    fHisto->SetHisto1D(2,fBinsP,p1,p2,1.0);
    fHisto->SetHisto1D(3,fBinsE,e1,e2,1.0);
    fHisto->SetHisto1D(4,fBinsE,e1,e2,1.0);
    fHisto->SetHisto1D(5,fBinsE,e1,e2,1.0);
    fHisto->SetHisto1D(6,fBinsE,e1,e2,1.0);
    fHisto->SetHisto1D(7,fBinsE,e1,e2,1.0);

  } else {
    fHisto->Add1D("h1","Elastic cross section (barn,1.0) as a functions of log10(p/GeV)",
                 fBinsP,p1,p2,1.0);
    fHisto->Add1D("h2","Elastic cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
    fHisto->Add1D("h3","Inelastic cross section (barn) as a functions of log10(p/GeV)",
                 fBinsP,p1,p2,1.0);
    fHisto->Add1D("h4","Inelastic cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
    fHisto->Add1D("h5","Capture cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
    fHisto->Add1D("h6","Fission cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
    fHisto->Add1D("h7","Charge exchange cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
    fHisto->Add1D("h8","Total cross section (barn) as a functions of log10(E/MeV)",
                 fBinsE,e1,e2,1.0);
  }

  fIsInitialised = true;
  fHisto->Book();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{
  if(fVerbose > 0) {
    G4cout << "HistoManager: End of run actions are started" << G4endl;
  }

  const G4Element* elm = 
    G4NistManager::Instance()->FindOrBuildElement(fElementName);
  const G4Material* mat = 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_"+fElementName);
  const G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);

  G4cout << "### Fill Cross Sections for " << fParticleName 
         << " off " << fElementName
         << G4endl;
  if(fVerbose > 0) {
    G4cout << "------------------------------------------------------------------------" 
           << G4endl;
    G4cout << "    N     E(MeV)   Elastic(b)   Inelastic(b)";
    if(particle == fNeutron) { G4cout << " Capture(b)   Fission(b)"; }
    G4cout << "   Total(b)" << G4endl;     
    G4cout << "------------------------------------------------------------------------" 
           << G4endl;
  }
  if(!particle || !elm) {
    G4cout << "HistoManager WARNING Particle or element undefined" << G4endl;
    return;
  }

  G4int prec = G4cout.precision();
  G4cout.precision(4);

  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  G4double mass = particle->GetPDGMass();

  // Build histograms

  G4double p1 = std::log10(fMinMomentum/GeV);
  G4double p2 = std::log10(fMaxMomentum/GeV);
  G4double e1 = std::log10(fMinKinEnergy/MeV);
  G4double e2 = std::log10(fMaxKinEnergy/MeV);
  G4double de = (e2 - e1)/G4double(fBinsE);
  G4double dp = (p2 - p1)/G4double(fBinsP);

  G4double x  = e1 - de*0.5; 
  G4double e, p, xs, xtot;
  G4int i;
  for(i=0; i<fBinsE; i++) {
    x += de;
    e  = std::pow(10.,x)*MeV;
    if(fVerbose>0) G4cout << std::setw(5) << i << std::setw(12) << e;  
    xs = store->GetElasticCrossSectionPerAtom(particle,e,elm,mat);
    xtot = xs;
    if(fVerbose>0) G4cout << std::setw(12) << xs/barn;  
    fHisto->Fill(1, x, xs/barn);    
    xs = store->GetInelasticCrossSectionPerAtom(particle,e,elm,mat);
    xtot += xs;
    if(fVerbose>0) G4cout << " " << std::setw(12) << xs/barn;  
    fHisto->Fill(3, x, xs/barn);    
    if(particle == fNeutron) {
      xs = store->GetCaptureCrossSectionPerAtom(particle,e,elm,mat);
      xtot += xs;
      if(fVerbose>0) G4cout << " " << std::setw(12) << xs/barn;  
      fHisto->Fill(4, x, xs/barn);    
      xs = store->GetFissionCrossSectionPerAtom(particle,e,elm,mat);
      xtot += xs;
      if(fVerbose>0) G4cout << " " << std::setw(12) << xs/barn;  
      fHisto->Fill(5, x, xs/barn);    
    }
    xs = store->GetChargeExchangeCrossSectionPerAtom(particle,e,elm,mat);
    if(fVerbose>0) G4cout << " " << std::setw(12) << xtot/barn << G4endl;   
    fHisto->Fill(6, x, xs/barn);    
    fHisto->Fill(7, x, xtot/barn);    
  }

  x = p1 - dp*0.5; 
  for(i=0; i<fBinsP; i++) {
    x += dp;
    p  = std::pow(10.,x)*GeV;
    e  = std::sqrt(p*p + mass*mass) - mass;
    xs = store->GetElasticCrossSectionPerAtom(particle,e,elm,mat);
    fHisto->Fill(0, x, xs/barn);    
    xs = store->GetInelasticCrossSectionPerAtom(particle,e,elm,mat);
    fHisto->Fill(2, x, xs/barn); 
  }
  if(fVerbose > 0) {
    G4cout << "-----------------------------------------------------------------" 
           << G4endl;
  }
  G4cout.precision(prec);
  fHisto->Save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  fVerbose = val; 
  fHisto->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

