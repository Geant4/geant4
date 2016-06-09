//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonGeneratorUniformSphere.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformSphere.cc,v 1.4 2006/06/29 16:16:35 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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
