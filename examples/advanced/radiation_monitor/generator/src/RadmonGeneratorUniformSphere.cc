//
// File name:     RadmonGeneratorUniformSphere.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformSphere.cc,v 1.1 2005-10-25 16:36:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorUniformSphere.hh"
#include "G4ParticleGun.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"

void                                            RadmonGeneratorUniformSphere :: ConvolveParticleGun(G4ParticleGun & gun)
{
 G4double radius(GetAttributeAsMeasure("Radius", "Length", -1.));
 
 if (radius<0)
 {
  G4cout << "RadmonGeneratorUniformSphere::ConvolveParticleGun: \"Radius\" not defined." << G4endl;
  return;
 }
 
 G4double phi(twopi*G4UniformRand());
 G4double cosTheta(2.*G4UniformRand()-1.);

 G4ThreeVector direction;
 direction.setRThetaPhi(1., std::acos(cosTheta), phi);
 
 G4RotationMatrix rotation;
 rotation.rotate(0., direction);
 
 direction*=radius;
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
 gun.SetParticlePosition(gun.GetParticlePosition()+direction);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorUniformSphere :: New(void) const
{
 return new RadmonGeneratorUniformSphere;
}
