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
// $Id: G4ParallelGeometrySampler.cc,v 1.1 2002-10-10 13:25:31 dressel Exp $
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
#include "G4VPScorer.hh"

#include "G4ParallelTransportConfigurator.hh"
#include "G4PScoreConfigurator.hh"
#include "G4PImportanceConfigurator.hh"
#include "G4WeightCutOffConfigurator.hh"

G4ParallelGeometrySampler::
G4ParallelGeometrySampler(G4VPhysicalVolume &worldvolume,
			  const G4String &particlename)
  :
  fParticleName(particlename),
  fParallelWorld(worldvolume),
  fParallelTransportConfigurator(0),
  fPImportanceConfigurator(0),
  fPScoreConfigurator(0),
  fWeightCutOffConfigurator(0),
  fIStore(0),
  fIsConfigured(false)
{}

G4ParallelGeometrySampler::~G4ParallelGeometrySampler(){
  ClearSampling();
};

void G4ParallelGeometrySampler::ClearSampling() {
  if (fParallelTransportConfigurator) delete fParallelTransportConfigurator;
  if (fPImportanceConfigurator) delete fPImportanceConfigurator;
  if (fPScoreConfigurator) delete fPScoreConfigurator;
  if (fWeightCutOffConfigurator) delete fWeightCutOffConfigurator;
  fIStore = 0;
};

G4bool G4ParallelGeometrySampler::IsConfigured() const{
  if (fIsConfigured) {
    G4cout << "G4ParallelGeometrySampler::CheckIfInit some initalization exists, use  ClearSampling() before a new initialization" << G4endl;
    return true;
  }
  return false;
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
}

void G4ParallelGeometrySampler::PrepareScoring(G4VPScorer *scorer){
  fPScoreConfigurator = 
    new G4PScoreConfigurator(fParticleName,
			     fParallelWorld.
			     GetParallelStepper(), 
			     *scorer);
}

void G4ParallelGeometrySampler::
PrepareWeightRoulett(G4double wsuvive, G4double wlimit, G4double isource){
  fWeightCutOffConfigurator = 
    new G4WeightCutOffConfigurator(fParticleName,
				   wsuvive,
				   wlimit,
				   isource,
				   fIStore,
				   &fParallelWorld.GetParallelStepper());
  
}

void G4ParallelGeometrySampler::Configure(){

  if (IsConfigured()) return;
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





