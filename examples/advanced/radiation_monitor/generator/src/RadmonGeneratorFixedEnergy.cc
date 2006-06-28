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
// File name:     RadmonGeneratorFixedEnergy.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedEnergy.cc,v 1.2 2006-06-28 13:53:45 gunter Exp $
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
