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
    void ConstructParticle() {};
  //
    void ConstructProcess(); 

  private:

  G4HadronElasticProcess  theElasticProcess;
  G4ProtonInelasticProcess theProtonInelastic;
  G4NeutronInelasticProcess  theNeutronInelastic;
  G4HadronElasticProcess* theNeutronElasticProcess;
  G4HadronFissionProcess* theFissionProcess;
  G4HadronCaptureProcess* theCaptureProcess;
  G4DeuteronInelasticProcess* theDeuteronInelasticProcess;
  G4TritonInelasticProcess* theTritonInelasticProcess;
  G4AlphaInelasticProcess* theAlphaInelasticProcess;
  G4IonInelasticProcess* theIonInelasticProcess;
};

#endif



