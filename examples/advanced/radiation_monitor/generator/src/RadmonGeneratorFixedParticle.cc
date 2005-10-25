//
// File name:     RadmonGeneratorFixedParticle.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedParticle.cc,v 1.1 2005-10-25 16:36:41 capra Exp $
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
