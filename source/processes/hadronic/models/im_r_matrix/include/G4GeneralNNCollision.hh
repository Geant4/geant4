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
#ifndef G4GeneralNNCollision_h
#define G4GeneralNNCollision_h

#include "G4CollisionComposite.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

class G4GeneralNNCollision : public G4CollisionComposite
{
  public:

  G4bool 
  IsInCharge(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
  {
    G4bool result = false;
    G4ParticleDefinition * aD = trk1.GetDefinition();
    G4ParticleDefinition * bD = trk2.GetDefinition();
    if(  (aD==G4Proton::Proton() || aD == G4Neutron::Neutron())
       &&(bD==G4Proton::Proton() || bD == G4Neutron::Neutron()) ) result = true;
    return result;
  }

};

#endif
