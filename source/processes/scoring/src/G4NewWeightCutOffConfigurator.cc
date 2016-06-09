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
// $Id: G4NewWeightCutOffConfigurator.cc,v 1.1 2006/11/20 10:02:16 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4NewWeightCutOffConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4NewWeightCutOffConfigurator.hh"
#include "G4NewWeightCutOffProcess.hh"

G4NewWeightCutOffConfigurator::
G4NewWeightCutOffConfigurator(G4VPhysicalVolume* worldvolume,
			      const G4String &particlename,
                                 G4double wsurvival,
                                 G4double wlimit,
                                 G4double isource,
                                 G4VIStore *istore,
                           const G4VGCellFinder &aGCellfinder, G4bool para)
  : fWorld(worldvolume),
    fPlacer(particlename),
    fPlaced(false),
    paraflag(para)
{
  fNewWeightCutOffProcess =
    new G4NewWeightCutOffProcess(wsurvival,wlimit,isource,istore,aGCellfinder,"NewWeightCutOffProcess",paraflag);
  if (!fNewWeightCutOffProcess)
  {
    G4Exception("G4NewWeightCutOffConfigurator::G4NewWeightCutOffConfigurator()",
                "FatalError", FatalException,
                "Failed to allocate G4NewWeightCutOffProcess !");
  }
}

G4NewWeightCutOffConfigurator::~G4NewWeightCutOffConfigurator()
{
  if (fPlaced)
  {
    fPlacer.RemoveProcess(fNewWeightCutOffProcess);
    delete fNewWeightCutOffProcess;
  }
}

void G4NewWeightCutOffConfigurator::Configure(G4VSamplerConfigurator *)
{
  G4cout << " entering new weight window configure " << G4endl;

  if(paraflag) fNewWeightCutOffProcess->SetParallelWorld(fWorld);

  fPlacer.AddProcessAsLastDoIt(fNewWeightCutOffProcess); 
  fPlaced = true;
}

const G4VTrackTerminator
*G4NewWeightCutOffConfigurator::GetTrackTerminator() const
{
  return 0;
}

