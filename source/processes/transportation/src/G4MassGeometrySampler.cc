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
// $Id: G4MassGeometrySampler.cc,v 1.1 2002-10-10 13:25:31 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassGeometrySampler.cc
//
// ----------------------------------------------------------------------



#include "G4MassGeometrySampler.hh"

#include "G4VIStore.hh"
#include "G4VPScorer.hh"

#include "G4MScoreConfigurator.hh"
#include "G4MImportanceConfigurator.hh"
#include "G4WeightCutOffConfigurator.hh"

G4MassGeometrySampler::
G4MassGeometrySampler(const G4String &particlename)
  :
  fParticleName(particlename),
  fMImportanceConfigurator(0),
  fMScoreConfigurator(0),
  fWeightCutOffConfigurator(0)
{}

G4MassGeometrySampler::~G4MassGeometrySampler(){
  ClearSampling();
};

void G4MassGeometrySampler::ClearSampling() {

  if (fMImportanceConfigurator) delete fMImportanceConfigurator;
  if (fMScoreConfigurator) delete fMScoreConfigurator;
  if (fWeightCutOffConfigurator) delete fWeightCutOffConfigurator;
  fIStore = 0;
  
};

G4bool G4MassGeometrySampler::IsConfigured() const{

  if (fIsConfigured) {
   G4cout << "G4MassGeometrySampler::CheckIfInit some initalization exists, use  ClearSampling() before a new initialization" << G4endl;
    return true;
  }
  return false;
}

void G4MassGeometrySampler::
PrepareScoring(G4VPScorer *scorer){
  fMScoreConfigurator = 
    new G4MScoreConfigurator(fParticleName, *scorer);
}

void G4MassGeometrySampler::
PrepareImportanceSampling(G4VIStore *istore,
			  const G4VImportanceAlgorithm 
			  *ialg) {
  fIStore = istore;

  fMImportanceConfigurator =
    new G4MImportanceConfigurator(fParticleName,
				    *fIStore,
				    ialg);
  
}

void G4MassGeometrySampler::
PrepareWeightRoulett(G4double wsurvive, 
		     G4double wlimit,
		     G4double isource) {
  fWeightCutOffConfigurator = 
    new G4WeightCutOffConfigurator(fParticleName,
				   wsurvive,
				   wlimit,
				   isource,
				   fIStore,
				   0);
  
}

void G4MassGeometrySampler::Configure(){
  if (IsConfigured()) return;
  fIsConfigured = true;

  if (fMScoreConfigurator) {
    fConfigurators.push_back(fMScoreConfigurator);
  }
  if (fMImportanceConfigurator) {
    fConfigurators.push_back(fMImportanceConfigurator);
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







