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

#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#include "G4InuclElementaryParticle.hh"
#endif

class G4ParticleLargerEkin {

public:
  
  G4bool operator() (const G4InuclElementaryParticle& part1,
		     const G4InuclElementaryParticle& part2) {

    return part1.getKineticEnergy() >= part2.getKineticEnergy();
    //  return part1.getEnergy() >= part2.getEnergy();
    //  return part1.getMomModule() >= part2.getMomModule();
  };
 
};

#endif // G4PARTICLE_LARGER_EKIN_HH
