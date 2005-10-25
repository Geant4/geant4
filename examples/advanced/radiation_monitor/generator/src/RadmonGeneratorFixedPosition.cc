//
// File name:     RadmonGeneratorFixedPosition.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedPosition.cc,v 1.1 2005-10-25 16:36:41 capra Exp $
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
