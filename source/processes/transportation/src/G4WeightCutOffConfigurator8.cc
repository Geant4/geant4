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
// $Id: G4WeightCutOffConfigurator8.cc,v 1.2 2007-06-01 07:53:27 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffConfigurator8
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4WeightCutOffConfigurator8.hh"
#include "G4WeightCutOffProcess8.hh"

G4WeightCutOffConfigurator8::
G4WeightCutOffConfigurator8(const G4String &particlename,
                                 G4double wsurvival,
                                 G4double wlimit,
                                 G4double isource,
                                 G4VIStore *istore,
                           const G4VGCellFinder &aGCellfinder)
  : fPlacer(particlename),
    fPlaced(false)
{
  fWeightCutOffProcess8 =
    new G4WeightCutOffProcess8(wsurvival,wlimit,isource,istore,aGCellfinder);
  if (!fWeightCutOffProcess8)
  {
    G4Exception("G4WeightCutOffConfigurator8::G4WeightCutOffConfigurator8()",
                "FatalError", FatalException,
                "Failed to allocate G4WeightCutOffProcess8 !");
  }
}

G4WeightCutOffConfigurator8::~G4WeightCutOffConfigurator8()
{
  if (fPlaced)
  {
    fPlacer.RemoveProcess(fWeightCutOffProcess8);
    delete fWeightCutOffProcess8;
  }
}

void G4WeightCutOffConfigurator8::Configure(G4VSamplerConfigurator8 *)
{
  fPlacer.AddProcessAsLastDoIt(fWeightCutOffProcess8); 
  fPlaced = true;
}

const G4VTrackTerminator
*G4WeightCutOffConfigurator8::GetTrackTerminator() const
{
  return 0;
}

