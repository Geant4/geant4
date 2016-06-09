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
//
// $Id: RE01BosonPhysics.cc,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//

#include "RE01BosonPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"


RE01BosonPhysics::RE01BosonPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{;}


RE01BosonPhysics::~RE01BosonPhysics()
{;}


void RE01BosonPhysics::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  // gamma
  G4Gamma::GammaDefinition();  
}


void RE01BosonPhysics::ConstructProcess()
{
   // Add e+e- pair creation, Compton scattering and photo-electric effect
   // to gamma

   G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(new G4GammaConversion());
   pManager->AddDiscreteProcess(new G4ComptonScattering());
   pManager->AddDiscreteProcess(new G4PhotoElectricEffect());
}


