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
// $Id: G4RangeTest.cc,v 1.6 2003-01-22 18:47:29 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 05 Oct 2001   MGP        Created
// 21 Jan 2003   VI         Cut per region
//
// -------------------------------------------------------------------

#include "G4RangeTest.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

G4RangeTest::~G4RangeTest()
{ }

G4bool G4RangeTest::Escape(const G4ParticleDefinition* particle,
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
    G4double rMin = G4std::min(cut,safety);
    value = (range > rMin);
  }
  return value;
}
