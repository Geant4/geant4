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
// $Id: G4ParallelGeometrySampler.cc,v 1.10 2003/11/26 14:51:50 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelGeometrySampler.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelGeometrySampler.hh"

#include "G4VIStore.hh"
#include "G4WeightWindowStore.hh"
#include "G4VScorer.hh"

#include "G4ParallelTransportConfigurator.hh"
#include "G4PScoreConfigurator.hh"
#include "G4PImportanceConfigurator.hh"
#include "G4WeightCutOffConfigurator.hh"
#include "G4ParallelGCellFinder.hh"
#include "G4PWeightWindowConfigurator.hh"

G4ParallelGeometrySampler::
G4ParallelGeometrySampler(G4VPhysicalVolume &worldvolume,
                          const G4String &particlename)
  : fParticleName(particlename),
    fParallelWorld(worldvolume),
    fParallelTransportConfigurator(0),
    fPImportanceConfigurator(0),
    fPScoreConfigurator(0),
    fGCellFinder(0),
    fWeightCutOffConfigurator(0),
    fIStore(0),
    fPWeightWindowConfigurator(0),
    fWWStore(0),
    fIsConfigured(false)
{
}

G4ParallelGeometrySampler::~G4ParallelGeometrySampler()
{
  ClearSampling();
}

void G4ParallelGeometrySampler::ClearSampling()
{
  if (fParallelTransportConfigurator)
  {
    delete fParallelTransportConfigurator;
    fParallelTransportConfigurator = 0;
  }
  if (fPWeightWindowConfigurator)
  {
    delete fPWeightWindowConfigurator;
    fPWeightWindowConfigurator = 0;
  }
  if (fPImportanceConfigurator)
  {
    delete fPImportanceConfigurator;
    fPImportanceConfigurator = 0;
  }
  if (fPScoreConfigurator)
  {
    delete fPScoreConfigurator;
    fPScoreConfigurator = 0;
  }
  if (fWeightCutOffConfigurator)
  {
    delete fWeightCutOffConfigurator;
    fWeightCutOffConfigurator = 0;
  }
  if (fGCellFinder)
  {
    delete fGCellFinder;
    fGCellFinder = 0;
  }
  fIStore = 0;
  fConfigurators.clear();
  fIsConfigured = false;
}

G4bool G4ParallelGeometrySampler::IsConfigured() const
{
  G4bool isconf = false;
  if (fIsConfigured)
  {
    G4cout << "WARNING - G4ParallelGeometrySampler::IsConfigured()"
           << "          Some initalization exists already."
           << "          Use ClearSampling() before a new initialization !"
           << G4endl;
    isconf =  true;
  }
  return isconf;
}

void G4ParallelGeometrySampler::
PrepareImportanceSampling(G4VIStore *istore,
                          const G4VImportanceAlgorithm *ialg)
{
  fIStore = istore;
  fPImportanceConfigurator = new G4PImportanceConfigurator(fParticleName,
                                                           fParallelWorld,
                                                           *istore, ialg);
  if (!fPImportanceConfigurator)
  {
    G4Exception("G4ParallelGeometrySampler::PrepareImportanceSampling()",
                "FatalError", FatalException,
                "Failed to create G4PImportanceConfigurator !");
  }
}

void G4ParallelGeometrySampler::PrepareScoring(G4VScorer *scorer)
{
  fPScoreConfigurator = 
    new G4PScoreConfigurator(fParticleName,
                             fParallelWorld.
                             GetParallelStepper(), 
                             *scorer);
  if (!fPScoreConfigurator)
  {
    G4Exception("G4ParallelGeometrySampler::PrepareScoring()",
                "FatalError", FatalException,
                "Failed to create G4PScoreConfigurator !");
  }
}

void
G4ParallelGeometrySampler::PrepareWeightRoulett(G4double wsuvive,
                                                G4double wlimit,
                                                G4double isource)
{
  fGCellFinder = new G4ParallelGCellFinder(fParallelWorld.
                                           GetParallelStepper());
  if (!fGCellFinder)
  {
    G4Exception("G4ParallelGeometrySampler::PrepareWeightRoulett()",
                "FatalError", FatalException,
                "Failed to allocate G4ParallelGCellFinder !");
  }

  fWeightCutOffConfigurator = 
    new G4WeightCutOffConfigurator(fParticleName,
                                   wsuvive,
                                   wlimit,
                                   isource,
                                   fIStore,
                                   *fGCellFinder);
  if (!fWeightCutOffConfigurator)
  {
    G4Exception("G4ParallelGeometrySampler::PrepareWeightRoulett()",
                "FatalError", FatalException,
                "Failed to allocate G4WeightCutOffConfigurator !");
  }
}

void
G4ParallelGeometrySampler::PrepareWeightWindow(G4VWeightWindowStore *wwstore,
                                               G4VWeightWindowAlgorithm *wwAlg,
                                               G4PlaceOfAction placeOfAction)
{
  if (placeOfAction != onBoundary)
  {
    // an additional paralel transport is needed
    // to cause a step on parallel boundaries
    fParallelTransportConfigurator = 
      new G4ParallelTransportConfigurator(fParticleName, fParallelWorld);    
  }
  
  fWWStore = wwstore;
  fPWeightWindowConfigurator = 
    new G4PWeightWindowConfigurator(fParticleName,
                                    fParallelWorld,
                                    *wwstore,
                                    wwAlg,
                                    placeOfAction);
}

void G4ParallelGeometrySampler::Configure()
{
  if (!IsConfigured())
  {
    fIsConfigured = true;

    if (fPScoreConfigurator)
    {
      fConfigurators.push_back(fPScoreConfigurator);
    }
    if (fPImportanceConfigurator)
    {
      fConfigurators.push_back(fPImportanceConfigurator);
    }
    else if (fPWeightWindowConfigurator)
    {
      fConfigurators.push_back(fPWeightWindowConfigurator);
      if (fParallelTransportConfigurator)
      {
        fConfigurators.push_back(fParallelTransportConfigurator);
      }
    }
    else
    {
      fParallelTransportConfigurator = 
        new G4ParallelTransportConfigurator(fParticleName, fParallelWorld);
      if (!fParallelTransportConfigurator)
      {
        G4Exception("G4ParallelGeometrySampler::Configure()",
                    "FatalError", FatalException,
                    "Failed to allocate G4ParallelTransportConfigurator !");
      }
      fConfigurators.push_back(fParallelTransportConfigurator);
    }
    G4VSamplerConfigurator *preConf = 0;
    for (G4Configurators::iterator it = fConfigurators.begin();
         it != fConfigurators.end(); it++)
    {
      G4VSamplerConfigurator *currConf =*it;
      currConf->Configure(preConf);
      preConf = *it;
    }
    
    if (fWeightCutOffConfigurator)
    {
      fWeightCutOffConfigurator->Configure(0);
    }
  }
  return;
}
