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
// Geant4 header G4EmBuilder
//
// Author V.Ivanchenko 22.05.2020
//
// EM physics for HEP hadrons and light ions
//

#ifndef G4EmBuilder_h
#define G4EmBuilder_h 1

#include "globals.hh"
#include <vector>

class G4hMultipleScattering;
class G4ParticleDefinition;
class G4NuclearStopping;
class G4VMscModel;

class G4EmBuilder
{
public:

  static void ConstructBasicEmPhysics(G4hMultipleScattering* hmsc, 
                                      const std::vector<G4int>& listHadrons);

  static void ConstructLightHadrons(G4ParticleDefinition* part1, 
                                    G4ParticleDefinition* part2,
                                    G4bool isHEP, G4bool isProton,
                                    G4bool isWVI);

  static void ConstructLightHadronsSS(G4ParticleDefinition* part1, 
                                      G4ParticleDefinition* part2,
				      G4bool isHEP);

  static void ConstructIonEmPhysics(G4hMultipleScattering* hmsc, 
                                    G4NuclearStopping* nucStopping); 

  static void ConstructIonEmPhysicsSS(); 

  // main method to be called from EM constructors to construct
  // EM physics for the list of leptons and hadrons common for
  // EM constructors
  static void ConstructCharged(G4hMultipleScattering* hmsc, 
                               G4NuclearStopping* nucStopping,
                               G4bool isWVI = true);

  static void ConstructChargedSS(G4hMultipleScattering* hmsc); 

  // minimal set of particles for EM physics
  static void ConstructMinimalEmSet();

  // prepare EM physics for construction
  static void PrepareEMPhysics();

  static void ConstructElectronMscProcess(G4VMscModel* msc1, G4VMscModel* msc2,
                                          G4ParticleDefinition* particle);
};

#endif


