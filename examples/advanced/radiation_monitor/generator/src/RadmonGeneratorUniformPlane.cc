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
// File name:     RadmonGeneratorUniformPlane.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformPlane.cc,v 1.5 2006/06/29 16:16:33 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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
