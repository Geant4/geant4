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
// $Id: A01EMPhysics.cc,v 1.7 2004/11/23 04:02:02 tkoi Exp $
// --------------------------------------------------------------
//
//
// 09-Oct-2003 Change gamma, electron, positorn process T. Koi
// 10-Jan-2004 Add Brems. of AlongStepDoIt for e+- T. Koi


#include "A01EMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>


A01EMPhysics::A01EMPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{
}

A01EMPhysics::~A01EMPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

#include "G4ProcessManager.hh"

void A01EMPhysics::ConstructProcess()
{
   G4ProcessManager * pManager = 0;

   //Gamma
   pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(new G4GammaConversion());
   pManager->AddDiscreteProcess(new G4ComptonScattering());
   pManager->AddDiscreteProcess(new G4PhotoElectricEffect());

   //Electorn
   pManager = G4Electron::Electron()->GetProcessManager();
   G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
   G4VProcess* theeminusIonisation         = new G4eIonisation();
   G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
   // 
   //  add process
   pManager->AddProcess(theeminusMultipleScattering);
   pManager->AddProcess(theeminusIonisation);
   pManager->AddProcess(theeminusBremsstrahlung);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
   pManager->SetProcessOrdering(theeminusBremsstrahlung,     idxAlongStep,3);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
   pManager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);

   //Positron
   pManager = G4Positron::Positron()->GetProcessManager();
   G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
   G4VProcess* theeplusIonisation         = new G4eIonisation();
   G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
   G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();

   pManager->AddProcess(theeplusMultipleScattering);
   pManager->AddProcess(theeplusIonisation);
   pManager->AddProcess(theeplusBremsstrahlung);
   pManager->AddProcess(theeplusAnnihilation);
   //
   // set ordering for AtRestDoIt
   pManager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
   //
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
   pManager->SetProcessOrdering(theeplusBremsstrahlung,     idxAlongStep,3);
   //
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
   pManager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
   pManager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);

}
