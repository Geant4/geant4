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
// File name:     RadmonGeneratorFixedParticle.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedParticle.cc,v 1.2 2006-06-28 13:53:47 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorFixedParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"

void                                            RadmonGeneratorFixedParticle :: ConvolveParticleGun(G4ParticleGun & gun)
{
 G4String particle(GetAttribute("Particle", "#"));
 
 if (particle=="#")
 {
  G4cout << "RadmonGeneratorFixedParticle::ConvolveParticleGun: \"Particle\" not defined." << G4endl;
  return;
 }

 G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
 G4ParticleDefinition * particleDefinition(particleTable->FindParticle(particle));
 
 if (particleDefinition==0)
 {
  G4cout << "RadmonGeneratorFixedParticle::ConvolveParticleGun: Particle \"" << particle << "\" not known." << G4endl;
  return;
 }
 
 gun.SetParticleDefinition(particleDefinition);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorFixedParticle :: New(void) const
{
 return new RadmonGeneratorFixedParticle;
}
