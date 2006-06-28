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
// File name:     RadmonGeneratorUniformSphere.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformSphere.cc,v 1.3 2006-06-28 13:54:01 gunter Exp $
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
 
 G4RotationMatrix rotation(G4RotationMatrix::IDENTITY);
 rotation.rotateY(direction.getTheta());
 rotation.rotateZ(phi);
 
 direction*=radius;
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
 gun.SetParticlePosition(gun.GetParticlePosition()+direction);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorUniformSphere :: New(void) const
{
 return new RadmonGeneratorUniformSphere;
}
