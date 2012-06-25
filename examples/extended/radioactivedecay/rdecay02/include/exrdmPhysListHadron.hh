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
/// \file radioactivedecay/rdecay02/include/exrdmPhysListHadron.hh
/// \brief Definition of the exrdmPhysListHadron class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef exrdmPhysListHadron_h
#define exrdmPhysListHadron_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4HadronElasticProcess.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4IonInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"


class exrdmPhysListHadron : public G4VPhysicsConstructor 
{
  public:
    exrdmPhysListHadron(const G4String& name = "hadron");
    virtual ~exrdmPhysListHadron();

  public:
  // Construct particle and physics
    virtual void ConstructParticle() {};
  //
    virtual void ConstructProcess();

  private:

  G4HadronElasticProcess  fTheElasticProcess;
  G4ProtonInelasticProcess fTheProtonInelastic;
  G4NeutronInelasticProcess  fTheNeutronInelastic;
  G4HadronElasticProcess* fTheNeutronElasticProcess;
  G4HadronFissionProcess* fTheFissionProcess;
  G4HadronCaptureProcess* fTheCaptureProcess;
  G4DeuteronInelasticProcess* fTheDeuteronInelasticProcess;
  G4TritonInelasticProcess* fTheTritonInelasticProcess;
  G4AlphaInelasticProcess* fTheAlphaInelasticProcess;
  G4IonInelasticProcess* fTheIonInelasticProcess;
};

#endif



