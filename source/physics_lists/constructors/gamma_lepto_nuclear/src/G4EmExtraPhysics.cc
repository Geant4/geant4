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
// $Id: G4EmExtraPhysics.cc 66704 2013-01-10 18:20:17Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 16.10.2012 A.Ribon: renamed G4EmExtraBertiniPhysics as G4EmExtraPhysics
// 10.04.2014 A.Dotti: Add MT functionality for messenger
// 24.04.2014 A.Ribon: switched on muon-nuclear by default
//
//----------------------------------------------------------------------------
//

#include "G4EmExtraPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4SynchrotronRadiation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4BertiniElectroNuclearBuilder.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4AutoDelete.hh"
 
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmExtraPhysics);

G4bool G4EmExtraPhysics::gnActivated  = true;
G4bool G4EmExtraPhysics::munActivated = true;
G4bool G4EmExtraPhysics::synActivated = false;
G4bool G4EmExtraPhysics::synActivatedForAll = false;

G4ThreadLocal G4BertiniElectroNuclearBuilder* G4EmExtraPhysics::theGNPhysics=0;
G4ThreadLocal G4SynchrotronRadiation* G4EmExtraPhysics::theSynchRad=0;

G4EmExtraPhysics::G4EmExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(ver)
{
  theMessenger = new G4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### G4EmExtraPhysics" << G4endl;
}

G4EmExtraPhysics::G4EmExtraPhysics(const G4String&): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(1)
{
  theMessenger = new G4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### G4EmExtraPhysics" << G4endl;
}

G4EmExtraPhysics::~G4EmExtraPhysics()
{
  delete theMessenger;
  theMessenger = 0;
}

void G4EmExtraPhysics::Synch(G4bool val)
{
  synActivated = val;
}

void G4EmExtraPhysics::SynchAll(G4bool val)
{
  synActivatedForAll = val;
  if(synActivatedForAll) { synActivated = true; }
}

void G4EmExtraPhysics::GammaNuclear(G4bool val)
{
  gnActivated = val;
}

void G4EmExtraPhysics::MuonNuclear(G4bool val)
{
  munActivated = val;
}

void G4EmExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

void G4EmExtraPhysics::ConstructProcess()
{
  G4ParticleDefinition* electron = G4Electron::Electron();
  G4ParticleDefinition* positron = G4Positron::Positron();
  G4ParticleDefinition* muonplus = G4MuonPlus::MuonPlus();
  G4ParticleDefinition* muonminus = G4MuonMinus::MuonMinus();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  if(gnActivated) {
    theGNPhysics = new G4BertiniElectroNuclearBuilder();
    theGNPhysics->Build();
    G4AutoDelete::Register(theGNPhysics);
  }
  if(munActivated) {
    G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
    G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
    muNucProcess->RegisterMe(muNucModel);
    ph->RegisterProcess( muNucProcess, muonplus);
    ph->RegisterProcess( muNucProcess, muonminus);
  }
  if(synActivated) {
    theSynchRad = new G4SynchrotronRadiation();
    ph->RegisterProcess( theSynchRad, electron);
    ph->RegisterProcess( theSynchRad, positron);
    G4AutoDelete::Register(theSynchRad);
    if(synActivatedForAll) {
      aParticleIterator->reset();
      G4ParticleDefinition* particle=0;

      while( (*aParticleIterator)() ) {
	particle = aParticleIterator->value();
	if( particle->GetPDGStable() && particle->GetPDGCharge() != 0.0) { 
	  if(verbose > 1) {
	    G4cout << "### G4SynchrotronRadiation for " 
		   << particle->GetParticleName() << G4endl;
	  }
	  ph->RegisterProcess( theSynchRad, particle);
	}
      }
    }
  }
}

