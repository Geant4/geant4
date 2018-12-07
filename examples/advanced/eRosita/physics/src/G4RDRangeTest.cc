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
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 05 Oct 2001   MGP        Created
// 21 Jan 2003   VI         Cut per region
//
// -------------------------------------------------------------------

#include "G4RDRangeTest.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

G4RDRangeTest::~G4RDRangeTest()
{ }

G4bool G4RDRangeTest::Escape(const G4ParticleDefinition* particle,
			   const G4MaterialCutsCouple* couple,
			   G4double energy,
			   G4double safety) const
{
  G4bool value = true;
  size_t idx = 0;
  if(particle == G4Electron::Electron()) idx = 1;
  else if(particle == G4Positron::Positron()) idx = 2;
  if(idx>0) {
    G4double range = G4EnergyLossTables::GetRange(particle,energy,couple);
    G4double cut = couple->GetProductionCuts()->GetProductionCut(idx);
    G4double rMin = std::min(cut,safety);
    value = (range > rMin);
  }
  return value;
}
