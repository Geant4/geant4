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
// $Id: G4PWeightWindowConfigurator.cc,v 1.4 2003/11/26 14:51:50 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// Class G4PWeightWindowConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4PWeightWindowConfigurator.hh"
#include "G4ParallelWWnTransportProcess.hh"
#include "G4ParallelWeightWindowProcess.hh" 
#include "G4WeightWindowAlgorithm.hh"
#include "G4ParallelWorld.hh"
#include "G4ParallelWeightWindowProcess.hh"

G4PWeightWindowConfigurator::
G4PWeightWindowConfigurator(const G4String &particlename,
                                  G4ParallelWorld &parallelWorld,
                                  G4VWeightWindowStore &wwstore,
                            const G4VWeightWindowAlgorithm *wwAlg,
                                  G4PlaceOfAction placeOfAction)
  : fPlacer(particlename),
    fPWorld(parallelWorld),
    fDeleteWWalg( ( ! wwAlg) ),
    fWWalgorithm(( (fDeleteWWalg) ? new G4WeightWindowAlgorithm(5,3,5)
                                  : wwAlg)),
    fExaminer(*fWWalgorithm, parallelWorld.GetParallelStepper(), wwstore),
    fParallelWWProcess(0),
    fPlaceOfAction(placeOfAction),
    fTrackTerminator(0)
{
}

G4PWeightWindowConfigurator::~G4PWeightWindowConfigurator()
{  
  if (fParallelWWProcess)
  {
    fPlacer.RemoveProcess(fParallelWWProcess);
    delete fParallelWWProcess;
    fTrackTerminator = 0;
    fParallelWWProcess = 0;
  }
  if (fDeleteWWalg)
  {
    delete fWWalgorithm;
  }
}

void
G4PWeightWindowConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  }

  if (fPlaceOfAction == onBoundary)
  {
    G4ParallelWWnTransportProcess *parallelWWnTransportProcess = 
      new G4ParallelWWnTransportProcess(fExaminer, 
                                        fPWorld.GetGeoDriver(), 
                                        fPWorld.GetParallelStepper(),
                                        terminator);
    fTrackTerminator = parallelWWnTransportProcess;
    fParallelWWProcess = parallelWWnTransportProcess;
  }
  else
  {
    G4ParallelWeightWindowProcess *parallelWeightWindowProcess = 
      new G4ParallelWeightWindowProcess(fExaminer, 
                                        fPWorld.GetParallelStepper(),
                                        terminator,
                                        fPlaceOfAction);
    fTrackTerminator = parallelWeightWindowProcess;
    fParallelWWProcess = parallelWeightWindowProcess;
  }

  fPlacer.AddProcessAsSecondDoIt(fParallelWWProcess);
}

const G4VTrackTerminator *
G4PWeightWindowConfigurator::GetTrackTerminator() const 
{
  return fTrackTerminator;
}
