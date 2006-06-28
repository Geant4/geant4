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
// File name:     RadmonGeneratorUniformPlane.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformPlane.cc,v 1.4 2006-06-28 13:53:59 gunter Exp $
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
 
 G4RotationMatrix rotation(GetAttributeAsRotationMatrix("Direction", G4RotationMatrix::IDENTITY));
 
 G4double x((G4UniformRand()-0.5)*width);
 G4double y((G4UniformRand()-0.5)*height);

 G4ThreeVector offset(x, y, 0.);
 
 offset.transform(rotation);
 
 gun.SetParticleMomentumDirection(rotation * gun.GetParticleMomentumDirection());
 gun.SetParticlePosition(gun.GetParticlePosition()+offset);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorUniformPlane :: New(void) const
{
 return new RadmonGeneratorUniformPlane;
}
