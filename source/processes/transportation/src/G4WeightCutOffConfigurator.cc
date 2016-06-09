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
// $Id: G4WeightCutOffConfigurator.cc,v 1.5 2003/11/26 14:51:50 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

