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
//
// $Id: G4WeightWindowConfigurator.cc 77938 2013-11-29 15:07:30Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4WeightWindowConfigurator.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4WeightWindowProcess.hh"

G4WeightWindowConfigurator::
G4WeightWindowConfigurator(const G4VPhysicalVolume* worldvolume,
			   const G4String &particlename,
                            G4VWeightWindowStore &wwstore,
                            const G4VWeightWindowAlgorithm *wwAlg,
                            G4PlaceOfAction placeOfAction, G4bool para)
  : fWorld(worldvolume),
    fPlacer(particlename),
    fWeightWindowStore(wwstore),
    fDeleteWWalg( ( ! wwAlg) ),
    fWWalgorithm(( (fDeleteWWalg) ? 
                   new G4WeightWindowAlgorithm(5,3,5) : wwAlg)),
    fWeightWindowProcess(0),
    fPlaceOfAction(placeOfAction)
{
  paraflag = para;
}

G4WeightWindowConfigurator::~G4WeightWindowConfigurator()
{  
  if (fWeightWindowProcess)
  {
    fPlacer.RemoveProcess(fWeightWindowProcess);
    delete fWeightWindowProcess;
  }
  if (fDeleteWWalg)
  {
    delete fWWalgorithm;
  }
}

void
G4WeightWindowConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  G4cout << " entering weight window configure " << G4endl;
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  };

  fWeightWindowProcess = 
    new G4WeightWindowProcess(*fWWalgorithm, 
                                  fWeightWindowStore, 
                                  terminator,
                                  fPlaceOfAction,"WeightWindowProcess",paraflag);

  if(paraflag) fWeightWindowProcess->SetParallelWorld(fWorld->GetName());

  fPlacer.AddProcessAsSecondDoIt(fWeightWindowProcess);
}

const G4VTrackTerminator *
G4WeightWindowConfigurator::GetTrackTerminator() const 
{
  return fWeightWindowProcess;
}

