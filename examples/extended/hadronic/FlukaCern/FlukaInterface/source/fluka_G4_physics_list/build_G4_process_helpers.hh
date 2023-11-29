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
// Helper: 
// Construct a G4HadronInelasticProcess (could be templated on process class),
// assign XS and model to the process, 
// and register the process to the process manager.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef BUILD_G4_PROCESS_HELPERS_HH
#define BUILD_G4_PROCESS_HELPERS_HH


#include <vector>

// G4
#include "globals.hh"


class G4ParticleDefinition;
class G4PhysicsListHelper;
class G4VCrossSectionDataSet;
class G4HadronicInteraction;


namespace build_G4_process_helpers {

  void buildInelasticProcessForEachParticle(const std::vector<G4int>& partList,
                                            G4PhysicsListHelper* ph,
                                            G4VCrossSectionDataSet* xs,
                                            //const G4double xsFactor,
                                            G4HadronicInteraction* hadronicModel);


  void buildInelasticProcess(G4ParticleDefinition* particle,
                             G4PhysicsListHelper* ph,
                             G4VCrossSectionDataSet* xs,
                             //const G4double xsFactor,
                             G4HadronicInteraction* hadronicModel);
}


#endif
#endif // G4_USE_FLUKA
