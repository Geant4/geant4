//
// File name:     RadmonGeneratorFixedDirection.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedDirection.cc,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorFixedDirection.hh"
#include "G4ParticleGun.hh"
#include "G4RotationMatrix.hh"

void                                            RadmonGeneratorFixedDirection :: ConvolveParticleGun(G4ParticleGun & gun)
{
 if (!ExistsAttribute("Direction"))
 {
  G4cout << "RadmonGeneratorFixedDirection::ConvolveParticleGun: \"Direction\" not defined." << G4endl;
  return;
 }
 
 G4ThreeVector direction(GetAttributeAsDirection("Directon", G4ThreeVector()));
 
 if (direction.getR()<0.5)
 {
  G4cout << "RadmonGeneratorFixedDirection::ConvolveParticleGun: \"Direction\" attribute has invalid format." << G4endl;
  return;
 }
 
 G4RotationMatrix rotation;
 rotation.rotate(0., direction);
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorFixedDirection :: New(void) const
{
 return new RadmonGeneratorFixedDirection;
}
