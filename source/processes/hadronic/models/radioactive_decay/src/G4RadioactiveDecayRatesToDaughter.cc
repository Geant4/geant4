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

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4RadioactiveDecayRatesToDaughter.hh"


G4RadioactiveDecayRatesToDaughter::G4RadioactiveDecayRatesToDaughter()
 : Z(0), A(0), E(0.0), generation(0), verboseLevel(0)
{}


G4RadioactiveDecayRatesToDaughter::
G4RadioactiveDecayRatesToDaughter(const G4RadioactiveDecayRatesToDaughter& right)
{
  Z = right.Z;
  A = right.A;
  E = right.E;
  generation = right.generation;
  decayRateC = right.decayRateC;
  taos = right.taos;
  verboseLevel = right.verboseLevel;
}


G4RadioactiveDecayRatesToDaughter& G4RadioactiveDecayRatesToDaughter::
operator=(const G4RadioactiveDecayRatesToDaughter& right)
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


G4RadioactiveDecayRatesToDaughter::~G4RadioactiveDecayRatesToDaughter()
{} 


void G4RadioactiveDecayRatesToDaughter::DumpInfo()
{
  G4cout << " Z: " << Z << "  A: " << A << "  E: " << E <<G4endl;
  G4cout << " Generation: " << generation << G4endl;
//  G4cout << " Coefficiency: " << decayRateC << endl;
//  G4cout << " Tao: " << tao << endl;
  // need to overload << for decayRAteC and tao first!

  G4cout << G4endl;
}







