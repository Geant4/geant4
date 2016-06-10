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
// $Id: G4WeightCutOffConfigurator.cc 77779 2013-11-28 07:49:08Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4WeightCutOffConfigurator.hh"
#include "G4WeightCutOffProcess.hh"

G4WeightCutOffConfigurator::
G4WeightCutOffConfigurator(const G4VPhysicalVolume* worldvolume,
			      const G4String &particlename,
                                 G4double wsurvival,
                                 G4double wlimit,
                                 G4double isource,
                                 G4VIStore *istore,
			         G4bool para)
  //                           const G4VGCellFinder &aGCellfinder, G4bool para)
  : fWorld(worldvolume),
    fPlacer(particlename),
    fPlaced(false),
    paraflag(para)
{
  fWeightCutOffProcess =
    new G4WeightCutOffProcess(wsurvival,wlimit,isource,istore,"WeightCutOffProcess",paraflag);
//     new G4WeightCutOffProcess(wsurvival,wlimit,isource,istore,aGCellfinder,"WeightCutOffProcess",paraflag);
  if (!fWeightCutOffProcess)
  {
    G4Exception("G4WeightCutOffConfigurator::G4WeightCutOffConfigurator()",
                "FatalError", FatalException,
                "Failed to allocate G4WeightCutOffProcess !");
  }
}

G4WeightCutOffConfigurator::~G4WeightCutOffConfigurator()
{
  if (fPlaced)
  {
    fPlacer.RemoveProcess(fWeightCutOffProcess);
    delete fWeightCutOffProcess;
  }
}

void G4WeightCutOffConfigurator::Configure(G4VSamplerConfigurator *)
{
  G4cout << " entering new weight window configure " << G4endl;

  if(paraflag) fWeightCutOffProcess->SetParallelWorld(fWorld->GetName());

  fPlacer.AddProcessAsLastDoIt(fWeightCutOffProcess); 
  fPlaced = true;
}

const G4VTrackTerminator
*G4WeightCutOffConfigurator::GetTrackTerminator() const
{
  return 0;
}

