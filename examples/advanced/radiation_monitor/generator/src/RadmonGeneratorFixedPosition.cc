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
// File name:     RadmonGeneratorFixedPosition.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedPosition.cc,v 1.2 2006-06-28 13:53:49 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorFixedPosition.hh"
#include "G4ParticleGun.hh"

void                                            RadmonGeneratorFixedPosition :: ConvolveParticleGun(G4ParticleGun & gun)
{
 if (!ExistsAttribute("Position"))
 {
  G4cout << "RadmonGeneratorFixedPosition::ConvolveParticleGun: \"Position\" not defined." << G4endl;
  return;
 }
 
 G4ThreeVector position(GetAttributeAsThreeVectorWithMeasure("Position", "Length", G4ThreeVector()));
 
 gun.SetParticlePosition(gun.GetParticlePosition()+position);
}



RadmonVGeneratorWithLabel *                     RadmonGeneratorFixedPosition :: New(void) const
{
 return new RadmonGeneratorFixedPosition;
}
