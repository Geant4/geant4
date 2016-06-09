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
// $Id: G4WeightCutOffConfigurator.cc,v 1.5.2.1 2006/06/29 21:12:48 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4WeightCutOffConfigurator.hh"
#include "G4WeightCutOffProcess.hh"

G4WeightCutOffConfigurator::
G4WeightCutOffConfigurator(const G4String &particlename,
                                 G4double wsurvival,
                                 G4double wlimit,
                                 G4double isource,
                                 G4VIStore *istore,
                           const G4VGCellFinder &aGCellfinder)
  : fPlacer(particlename),
    fPlaced(false)
{
  fWeightCutOffProcess =
    new G4WeightCutOffProcess(wsurvival,wlimit,isource,istore,aGCellfinder);
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
  fPlacer.AddProcessAsLastDoIt(fWeightCutOffProcess); 
  fPlaced = true;
}

const G4VTrackTerminator
*G4WeightCutOffConfigurator::GetTrackTerminator() const
{
  return 0;
}

