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
// $Id: HistoManager.cc,v 1.7 2010-10-11 11:02:36 vnivanch Exp $
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
  histo     = new Histo(verbose);
  neutron   = G4Neutron::Neutron();
  verbose   = 1;

  particleName  = "proton";
  elementName   = "Al";

  minKinEnergy  = 0.1*MeV;
  maxKinEnergy  = 10*TeV;
  minMomentum   = 1*MeV;
  maxMomentum   = 10*TeV;

  nBinsE    = 800;
  nBinsP    = 700;

  isInitialised = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

HistoManager::~HistoManager()
{
  delete histo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::BeginOfRun()
{
  G4double p1 = std::log10(minMomentum/GeV);
  G4double p2 = std::log10(maxMomentum/GeV);
  G4double e1 = std::log10(minKinEnergy/MeV);
  G4double e2 = std::log10(maxKinEnergy/MeV);

  //G4cout<<"e1= "<<e1<<" e2= "<<e2<<" p1= "<<p1<<" p2= "<<p2<<G4endl;

  if(isInitialised) {
    histo->setHisto1D(0,nBinsP,p1,p2);
    histo->setHisto1D(1,nBinsE,e1,e2);
    histo->setHisto1D(2,nBinsP,p1,p2);
    histo->setHisto1D(3,nBinsE,e1,e2);
    histo->setHisto1D(4,nBinsE,e1,e2);
    histo->setHisto1D(5,nBinsE,e1,e2);
    histo->setHisto1D(6,nBinsE,e1,e2);
    histo->setHisto1D(7,nBinsE,e1,e2);

  } else {
    histo->add1D("h1","Elastic cross section (barn) as a functions of log10(p/GeV)",
		 nBinsP,p1,p2);
    histo->add1D("h2","Elastic cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
    histo->add1D("h3","Inelastic cross section (barn) as a functions of log10(p/GeV)",
		 nBinsP,p1,p2);
    histo->add1D("h4","Inelastic cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
    histo->add1D("h5","Capture cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
    histo->add1D("h6","Fission cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
    histo->add1D("h7","Charge exchange cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
    histo->add1D("h8","Total cross section (barn) as a functions of log10(E/MeV)",
		 nBinsE,e1,e2);
  }

  isInitialised = true;
  histo->book();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::EndOfRun()
{
  if(verbose > 0) {
    G4cout << "HistoManager: End of run actions are started" << G4endl;
  }

  const G4Element* elm = 
    G4NistManager::Instance()->FindOrBuildElement(elementName);
  const G4ParticleDefinition* particle = 
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);

  G4cout << "### Fill Cross Sections for " << particleName 
	 << " off " << elementName
	 << G4endl;
  if(verbose > 0) {
    G4cout << "------------------------------------------------------------------------" 
	   << G4endl;
    G4cout << "    N     E(MeV)   Elastic(b)   Inelastic(b)";
    if(particle == neutron) { G4cout << " Capture(b)   Fission(b)"; }
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

  G4double p1 = std::log10(minMomentum/GeV);
  G4double p2 = std::log10(maxMomentum/GeV);
  G4double e1 = std::log10(minKinEnergy/MeV);
  G4double e2 = std::log10(maxKinEnergy/MeV);
  G4double de = (e2 - e1)/G4double(nBinsE);
  G4double dp = (p2 - p1)/G4double(nBinsP);

  G4double x  = e1 - de*0.5; 
  G4double e, p, xs, xtot;
  G4int i;
  for(i=0; i<nBinsE; i++) {
    x += de;
    e  = std::pow(10.,x)*MeV;
    if(verbose>0) G4cout << std::setw(5) << i << std::setw(12) << e;  
    xs = store->GetElasticCrossSectionPerAtom(particle,e,elm);
    xtot = xs;
    if(verbose>0) G4cout << std::setw(12) << xs/barn;  
    histo->fill(1, x, xs/barn);    
    xs = store->GetInelasticCrossSectionPerAtom(particle,e,elm);
    xtot += xs;
    if(verbose>0) G4cout << " " << std::setw(12) << xs/barn;  
    histo->fill(3, x, xs/barn);    
    if(particle == neutron) {
      xs = store->GetCaptureCrossSectionPerAtom(particle,e,elm);
      xtot += xs;
      if(verbose>0) G4cout << " " << std::setw(12) << xs/barn;  
      histo->fill(4, x, xs/barn);    
      xs = store->GetFissionCrossSectionPerAtom(particle,e,elm);
      xtot += xs;
      if(verbose>0) G4cout << " " << std::setw(12) << xs/barn;  
      histo->fill(5, x, xs/barn);    
    }
    xs = store->GetChargeExchangeCrossSectionPerAtom(particle,e,elm);
    if(verbose>0) G4cout << " " << std::setw(12) << xtot/barn << G4endl;   
    histo->fill(6, x, xs/barn);    
    histo->fill(7, x, xtot/barn);    
  }

  x = p1 - dp*0.5; 
  for(i=0; i<nBinsP; i++) {
    x += dp;
    p  = std::pow(10.,x)*GeV;
    e  = std::sqrt(p*p + mass*mass) - mass;
    xs = store->GetElasticCrossSectionPerAtom(particle,e,elm);
    histo->fill(0, x, xs/barn);    
    xs = store->GetInelasticCrossSectionPerAtom(particle,e,elm);
    histo->fill(2, x, xs/barn); 
  }
  if(verbose > 0) {
    G4cout << "-----------------------------------------------------------------" 
	   << G4endl;
  }
  G4cout.precision(prec);
  if(verbose > 1) { histo->print(0); }
  histo->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void HistoManager::SetVerbose(G4int val)        
{
  verbose = val; 
  histo->setVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

