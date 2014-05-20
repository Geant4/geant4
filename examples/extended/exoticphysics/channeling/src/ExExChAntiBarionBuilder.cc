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

#include "ExExChAntiBarionBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Wrapper
#include "XWrapperDiscreteProcess.hh"

ExExChAntiBarionBuilder::
ExExChAntiBarionBuilder(): wasActivated(false)
{
    theAntiProtonInelastic=new G4AntiProtonInelasticProcess;
    theAntiNeutronInelastic=new G4AntiNeutronInelasticProcess;
    theAntiDeuteronInelastic=new G4AntiDeuteronInelasticProcess;
    theAntiTritonInelastic=new G4AntiTritonInelasticProcess;
    theAntiHe3Inelastic=new G4AntiHe3InelasticProcess;
    theAntiAlphaInelastic=new G4AntiAlphaInelasticProcess;
}

ExExChAntiBarionBuilder::
~ExExChAntiBarionBuilder(){
    delete theAntiProtonInelastic;
    delete theAntiNeutronInelastic;
    delete theAntiDeuteronInelastic;
    delete theAntiTritonInelastic;
    delete theAntiHe3Inelastic;
    delete theAntiAlphaInelastic;
}

void ExExChAntiBarionBuilder::
Build()
{
    wasActivated = true;
    
    std::vector<G4VAntiBarionBuilder *>::iterator i;
    for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
    {
        (*i)->Build(theAntiProtonInelastic);
        (*i)->Build(theAntiNeutronInelastic);
        (*i)->Build(theAntiDeuteronInelastic);
        (*i)->Build(theAntiTritonInelastic);
        (*i)->Build(theAntiHe3Inelastic);
        (*i)->Build(theAntiAlphaInelastic);
    }
    G4ProcessManager * theProcMan;
    theProcMan = G4AntiProton::AntiProton()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiProtonInelastic_wr =
        new XWrapperDiscreteProcess();
    theAntiProtonInelastic_wr->RegisterProcess(theAntiProtonInelastic,1);
    theProcMan->AddDiscreteProcess(theAntiProtonInelastic_wr);
    
    theProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
    theProcMan->AddDiscreteProcess(theAntiNeutronInelastic);
    
    theProcMan = G4AntiDeuteron::AntiDeuteron()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiDeuteronInelastic_wr =
        new XWrapperDiscreteProcess();
    theAntiDeuteronInelastic_wr->RegisterProcess(theAntiDeuteronInelastic,1);
    theProcMan->AddDiscreteProcess(theAntiDeuteronInelastic_wr);
    
    theProcMan = G4AntiTriton::AntiTriton()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiTritonInelastic_wr =
        new XWrapperDiscreteProcess();
    theAntiTritonInelastic_wr->RegisterProcess(theAntiTritonInelastic,1);
    theProcMan->AddDiscreteProcess(theAntiTritonInelastic_wr);
    
    theProcMan = G4AntiHe3::AntiHe3()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiHe3Inelastic_wr =
        new XWrapperDiscreteProcess();
    theAntiHe3Inelastic_wr->RegisterProcess(theAntiHe3Inelastic,1);
    theProcMan->AddDiscreteProcess(theAntiHe3Inelastic_wr);
    
    theProcMan = G4AntiAlpha::AntiAlpha()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiAlphaInelastic_wr =
        new XWrapperDiscreteProcess();
    theAntiAlphaInelastic_wr->RegisterProcess(theAntiAlphaInelastic,1);
    theProcMan->AddDiscreteProcess(theAntiAlphaInelastic_wr);
}

