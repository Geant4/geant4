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
#ifdef G4_USE_FLUKA


#include "build_G4_process_helpers.hh"

// G4
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"
//#include "G4HadronicParameters.hh"
//#include "G4ProcessManager.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronInelasticProcess.hh"


// ***************************************************************************
// Helpers: 
// Construct a G4HadronInelasticProcess (could be templated on process class),
// assign XS and model to the process, 
// and register the process to the process manager.
// ***************************************************************************
namespace build_G4_process_helpers {

  // ***************************************************************************
  // For all particles.
  // ***************************************************************************
  void buildInelasticProcessForEachParticle(const std::vector<G4int>& allParticlesPDGIds,
                                            G4PhysicsListHelper* const helper,
                                            G4VCrossSectionDataSet* const xs,
                                            /*const G4double xsFactor,*/
                                            G4HadronicInteraction* const hadronicModel) {

    const auto particlesTable = G4ParticleTable::GetParticleTable();

    // Loop on all particles
    for (const auto& particlePDGId : allParticlesPDGIds) {

      const auto particle = particlesTable->FindParticle(particlePDGId);
      if (!particle) { continue; }

      buildInelasticProcess(particle,
                            helper,
                            xs,
                            //xsFactor,
                            hadronicModel);
    }
  }


  // ***************************************************************************
  // For a specific particle.
  // ***************************************************************************
  void buildInelasticProcess(G4ParticleDefinition* const particle,
                             G4PhysicsListHelper* const helper,
                             G4VCrossSectionDataSet* const xs,
                             //const G4double xsFactor,
                             G4HadronicInteraction* const hadronicModel) {

    const auto process = new G4HadronInelasticProcess(particle->GetParticleName() + "HadronInelastic", particle);
    process->AddDataSet(xs);
    process->RegisterMe(hadronicModel);
    //if (param->ApplyFactorXS()) ((G4HadronicProcess*)process)->MultiplyCrossSectionBy(xsFactor);
    helper->RegisterProcess(process, particle);
  }

} // namespace build_G4_process_helpers


#endif // G4_USE_FLUKA
