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
// $Id: G4RangeTest.cc,v 1.5 2002-05-28 09:20:21 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 05 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4RangeTest.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"

G4RangeTest::~G4RangeTest()
{ }

G4bool G4RangeTest::Escape(const G4ParticleDefinition* particle, 
			   const G4Material* material,
			   G4double energy, 
			   G4double safety) const
{
  G4double range = G4EnergyLossTables::GetRange(particle,energy,material);
  G4double cut = particle->GetRangeThreshold(material);
  G4double rMin = G4std::min(cut,safety);
  G4bool value = (range > rMin);

  return value;
}
