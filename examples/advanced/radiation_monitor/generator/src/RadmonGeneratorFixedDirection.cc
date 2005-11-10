//
// File name:     RadmonGeneratorFixedDirection.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedDirection.cc,v 1.2 2005-11-10 08:11:26 capra Exp $
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
 
 G4RotationMatrix rotation(GetAttributeAsRotationMatrix("Direction", G4RotationMatrix::IDENTITY));
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorFixedDirection :: New(void) const
{
 return new RadmonGeneratorFixedDirection;
}
