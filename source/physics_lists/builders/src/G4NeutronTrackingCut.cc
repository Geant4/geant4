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
// $Id: G4NeutronTrackingCut.cc,v 1.4 2006/11/24 16:30:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronTrackingCut
//
// Author: Nov 2006 G.Folger
//
//
//----------------------------------------------------------------------------
//

#include "G4NeutronTrackingCut.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4Neutron.hh"
#include "G4NeutronKiller.hh"

G4NeutronTrackingCut::G4NeutronTrackingCut(const G4String& name, G4int ver)
  :  G4VPhysicsConstructor(name), verbose(ver), wasActivated(false)
{}

G4NeutronTrackingCut::~G4NeutronTrackingCut()
{
  if(wasActivated) 
  {
      delete pNeutronKiller;
  }    
}

void G4NeutronTrackingCut::ConstructParticle()
{

// G4cout << "G4NeutronTrackingCut::ConstructParticle" << G4endl;

  G4Neutron::NeutronDefinition();
}

void G4NeutronTrackingCut::ConstructProcess()
{
  G4double timeLimit     =10*microsecond;
//  G4double kineticEnergylimit =1*keV;
  if(wasActivated) return;
  wasActivated = true;

  // Add Decay Process
  pNeutronKiller = new G4NeutronKiller();
  G4ParticleDefinition * particle = G4Neutron::Neutron();
  G4ProcessManager * pmanager = particle->GetProcessManager();
//      if(verbose > 1)
        G4cout << "### Adding Neutron tracking cut for " << particle->GetParticleName() << G4endl;
	G4cout << "###  cut value is " << timeLimit/microsecond << " microseconds" <<  G4endl;
      pmanager -> AddDiscreteProcess(pNeutronKiller);
//      pNeutronKiller->SetKinEnergyLimit(kineticEnergylimit);
      pNeutronKiller->SetTimeLimit(timeLimit);
}


