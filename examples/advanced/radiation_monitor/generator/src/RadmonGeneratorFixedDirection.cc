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
// File name:     RadmonGeneratorFixedDirection.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedDirection.cc,v 1.3 2006-06-28 13:53:43 gunter Exp $
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
