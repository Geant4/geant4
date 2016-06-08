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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// G4RadioactiveDecayRate.cc

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4RadioactiveDecayRate.hh"


G4RadioactiveDecayRate::G4RadioactiveDecayRate()
{
  ;
  //do nothing at the momment
}



G4RadioactiveDecayRate::G4RadioactiveDecayRate(const G4RadioactiveDecayRate &right)
{
  Z = right.Z;
  A = right.A;
  E = right.E;
  generation = right.generation;
  decayRateC = right.decayRateC;
  taos = right.taos;
  //  verboseLevel = right.verboseLevel;
}

G4RadioactiveDecayRate & G4RadioactiveDecayRate::operator=(const G4RadioactiveDecayRate &right)
{
  if (this != &right) { 
    Z = right.Z;
    A = right.A;
    E = right.E;
    generation = right.generation;
    decayRateC = right.decayRateC;
    taos = right.taos;
    //    verboseLevel = right.verboseLevel;
  }
  return *this;
}


G4RadioactiveDecayRate::~G4RadioactiveDecayRate()
{ ;} 


void G4RadioactiveDecayRate::DumpInfo()
{
  G4cout << " Z: " << Z << "  A: " << A << "  E: " << E <<G4endl;
  G4cout << " Generation: " << generation << G4endl;
//  G4cout << " Coefficiency: " << decayRateC << endl;
//  G4cout << " Tao: " << tao << endl;
  // need to overload << for decayRAteC and tao first!

  G4cout << G4endl;
}







