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

#include "ExExChPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Wrapper
#include "XWrapperDiscreteProcess.hh"

ExExChPiKBuilder::
ExExChPiKBuilder(): wasActivated(false)
{
    thePionPlusInelastic=new G4PionPlusInelasticProcess;
    thePionMinusInelastic=new G4PionMinusInelasticProcess;
    theKaonPlusInelastic=new G4KaonPlusInelasticProcess;
    theKaonMinusInelastic=new G4KaonMinusInelasticProcess;
    theKaonZeroLInelastic=new G4KaonZeroLInelasticProcess;
    theKaonZeroSInelastic=new G4KaonZeroSInelasticProcess;
}

ExExChPiKBuilder::
~ExExChPiKBuilder(){
    delete thePionPlusInelastic;
    delete thePionMinusInelastic;
    delete theKaonPlusInelastic;
    delete theKaonMinusInelastic;
    delete theKaonZeroLInelastic;
    delete theKaonZeroSInelastic;
}

void ExExChPiKBuilder::
Build()
{
    wasActivated = true;
    
    std::vector<G4VPiKBuilder *>::iterator i;
    for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
    {
        (*i)->Build(thePionPlusInelastic);
        (*i)->Build(thePionMinusInelastic);
        (*i)->Build(theKaonPlusInelastic);
        (*i)->Build(theKaonMinusInelastic);
        (*i)->Build(theKaonZeroLInelastic);
        (*i)->Build(theKaonZeroSInelastic);
    }
    G4ProcessManager * theProcMan;
    theProcMan = G4PionPlus::PionPlus()->GetProcessManager();
    XWrapperDiscreteProcess* thePionPlusInelastic_wrapper =
        new XWrapperDiscreteProcess();
    thePionPlusInelastic_wrapper->RegisterProcess(thePionPlusInelastic,1);
    theProcMan->AddDiscreteProcess(thePionPlusInelastic_wrapper);
    
    theProcMan = G4PionMinus::PionMinus()->GetProcessManager();
    XWrapperDiscreteProcess* thePionMinusInelastic_wrapper =
        new XWrapperDiscreteProcess();
    thePionMinusInelastic_wrapper->RegisterProcess(thePionMinusInelastic,1);
    theProcMan->AddDiscreteProcess(thePionMinusInelastic_wrapper);
    
    theProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
    XWrapperDiscreteProcess* theKaonPlusInelastic_wrapper =
        new XWrapperDiscreteProcess();
    theKaonPlusInelastic_wrapper->RegisterProcess(theKaonPlusInelastic,1);
    theProcMan->AddDiscreteProcess(theKaonPlusInelastic_wrapper);
    
    theProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theKaonMinusInelastic_wrapper =
        new XWrapperDiscreteProcess();
    theKaonMinusInelastic_wrapper->RegisterProcess(theKaonMinusInelastic,1);
    theProcMan->AddDiscreteProcess(theKaonMinusInelastic_wrapper);
    
    theProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
    theProcMan->AddDiscreteProcess(theKaonZeroLInelastic);
    
    theProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
    theProcMan->AddDiscreteProcess(theKaonZeroSInelastic);
}
// 2002 by J.P. Wellisch
