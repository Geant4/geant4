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
// $Id: G4ParallelGeometrySampler.cc,v 1.6 2002-11-04 10:47:56 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelGeometrySampler.cc
//
// ----------------------------------------------------------------------


#include "G4ParallelGeometrySampler.hh"

#include "G4VIStore.hh"
#include "G4VScorer.hh"

#include "G4ParallelTransportConfigurator.hh"
#include "G4PScoreConfigurator.hh"
#include "G4PImportanceConfigurator.hh"
#include "G4WeightCutOffConfigurator.hh"
#include "G4ParallelGCellFinder.hh"

G4ParallelGeometrySampler::
G4ParallelGeometrySampler(G4VPhysicalVolume &worldvolume,
			  const G4String &particlename)
  :
  fParticleName(particlename),
  fParallelWorld(worldvolume),
  fParallelTransportConfigurator(0),
  fPImportanceConfigurator(0),
  fPScoreConfigurator(0),
  fGCellFinder(0),
  fWeightCutOffConfigurator(0),
  fIStore(0),
  fIsConfigured(false)
{}

G4ParallelGeometrySampler::~G4ParallelGeometrySampler(){
  ClearSampling();
};

void G4ParallelGeometrySampler::ClearSampling() {
  if (fParallelTransportConfigurator) {
    delete fParallelTransportConfigurator;
    fParallelTransportConfigurator = 0;
  }
  if (fPImportanceConfigurator) {
    delete fPImportanceConfigurator;
    fPImportanceConfigurator = 0;
  }
  if (fPScoreConfigurator) {
    delete fPScoreConfigurator;
    fPScoreConfigurator = 0;
  }
  if (fWeightCutOffConfigurator) {
    delete fWeightCutOffConfigurator;
    fWeightCutOffConfigurator = 0;
  }
  if (fGCellFinder) {
    delete fGCellFinder;
    fGCellFinder = 0;
  }
  fIStore = 0;
  fConfigurators.clear();
  fIsConfigured = false;
};

G4bool G4ParallelGeometrySampler::IsConfigured() const{
  G4bool isconf = false;
  if (fIsConfigured) {
    G4cout << "G4ParallelGeometrySampler::CheckIfInit some initalization exists, use  ClearSampling() before a new initialization" << G4endl;
    isconf =  true;
  }
  return isconf;
}

void G4ParallelGeometrySampler::
PrepareImportanceSampling(G4VIStore *istore,
			  const G4VImportanceAlgorithm 
			  *ialg) {
  fIStore = istore;
  fPImportanceConfigurator = 
    new G4PImportanceConfigurator(fParticleName,
				  fParallelWorld,
				  *istore,
				  ialg);
  if (!fPImportanceConfigurator) {
    G4Exception("ERROR:G4ParallelGeometrySampler::PrepareImportanceSampling: new failed to create G4PImportanceConfigurator!");
  }
}

void G4ParallelGeometrySampler::PrepareScoring(G4VScorer *scorer){
  fPScoreConfigurator = 
    new G4PScoreConfigurator(fParticleName,
			     fParallelWorld.
			     GetParallelStepper(), 
			     *scorer);
  if (!fPScoreConfigurator) {
    G4Exception("ERROR:G4ParallelGeometrySampler::PrepareScoring: new failed to create G4PScoreConfigurator!");
  }
  
}

void G4ParallelGeometrySampler::
PrepareWeightRoulett(G4double wsuvive, G4double wlimit, G4double isource){
  fGCellFinder = new G4ParallelGCellFinder(fParallelWorld.
					   GetParallelStepper());
  if (!fGCellFinder) {
    G4Exception("ERROR:G4ParallelGeometrySampler::PrepareWeightRoulett: new failed to create G4ParallelGCellFinder!");
  }

  fWeightCutOffConfigurator = 
    new G4WeightCutOffConfigurator(fParticleName,
				   wsuvive,
				   wlimit,
				   isource,
				   fIStore,
				   *fGCellFinder);
  if (!fWeightCutOffConfigurator) {
    G4Exception("ERROR:G4ParallelGeometrySampler::PrepareWeightRoulett: new failed to create G4WeightCutOffConfigurator!");
  }
  
}

void G4ParallelGeometrySampler::Configure(){

  if (!IsConfigured()) {
    fIsConfigured = true;

    if (fPScoreConfigurator) {
      fConfigurators.push_back(fPScoreConfigurator);
    }
    if (fPImportanceConfigurator) {
      fConfigurators.push_back(fPImportanceConfigurator);
    }
    else {
      fParallelTransportConfigurator = 
	new G4ParallelTransportConfigurator(fParticleName,
					    fParallelWorld);
      if (!fParallelTransportConfigurator) {
	G4Exception("ERROR:G4ParallelGeometrySampler::Configure: new failed to create G4ParallelTransportConfigurator!");
      }
      fConfigurators.push_back(fParallelTransportConfigurator);
    }
    G4VSamplerConfigurator *preConf = 0;
    for (G4Configurators::iterator it = fConfigurators.begin();
	 it != fConfigurators.end(); it++) {
      G4VSamplerConfigurator *currConf =*it;
      currConf->Configure(preConf);
      preConf = *it;
    }
    
    if (fWeightCutOffConfigurator) {
      fWeightCutOffConfigurator->Configure(0);
    }
  }
  return;
}





