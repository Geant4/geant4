//
// File name:     RadmonGeneratorUniformPlane.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformPlane.cc,v 1.2 2005-11-10 08:11:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorUniformPlane.hh"
#include "G4ParticleGun.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"

void                                            RadmonGeneratorUniformPlane :: ConvolveParticleGun(G4ParticleGun & gun)
{
 G4double width(GetAttributeAsMeasure("Width", "Length", -1.)); // x
 if (width<0)
 {
  G4cout << "RadmonGeneratorUniformPlane::ConvolveParticleGun: \"Width\" not defined." << G4endl;
  return;
 }
 
 G4double height(GetAttributeAsMeasure("Height", "Length", -1.)); // y
 if (height<0)
 {
  G4cout << "RadmonGeneratorUniformPlane::ConvolveParticleGun: \"Height\" not defined." << G4endl;
  return;
 }
 
 if (!ExistsAttribute("Direction"))
 {
  G4cout << "RadmonGeneratorUniformPlane::ConvolveParticleGun: \"Direction\" not defined." << G4endl;
  return;
 }
 
 G4ThreeVector direction(GetAttributeAsDirection("Directon", G4ThreeVector()));
 
 if (direction.getR()<0.5)
 {
  G4cout << "RadmonGeneratorUniformPlane::ConvolveParticleGun: \"Direction\" attribute has invalid format." << G4endl;
  return;
 }
 
 G4double delta(GetAttributeAsMeasure("Delta", "Angle", -1.));
 if (delta<0)
 {
  G4cout << "RadmonGeneratorUniformPlane::ConvolveParticleGun: \"Delta\" not defined." << G4endl;
  return;
 }
 
 G4double x((G4UniformRand()-0.5)*width);
 G4double y((G4UniformRand()-0.5)*height);

 G4ThreeVector offset(x, y, 0.);
 
 G4RotationMatrix rotation;
 rotation.rotate(delta/rad, direction);
 
 offset.transform(rotation);
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
 gun.SetParticlePosition(gun.GetParticlePosition()+offset);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorUniformPlane :: New(void) const
{
 return new RadmonGeneratorUniformPlane;
}
