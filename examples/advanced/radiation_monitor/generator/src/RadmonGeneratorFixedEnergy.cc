//
// File name:     RadmonGeneratorFixedEnergy.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedEnergy.cc,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorFixedEnergy.hh"
#include "G4ParticleGun.hh"

void                                            RadmonGeneratorFixedEnergy :: ConvolveParticleGun(G4ParticleGun & gun)
{
 G4double energy(GetAttributeAsMeasure("Energy", "Energy", -1.));
 
 if (energy<0)
 {
  G4cout << "RadmonGeneratorFixedEnergy::ConvolveParticleGun: \"Energy\" not defined." << G4endl;
  return;
 }
 
 gun.SetParticleEnergy(gun.GetParticleEnergy()+energy);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorFixedEnergy :: New(void) const
{
 return new RadmonGeneratorFixedEnergy;
}
