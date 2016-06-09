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
// $Id: TiaraPrimaryGeneratorAction.cc,v 1.5 2006/06/29 15:45:23 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "globals.hh"

#include "TiaraPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Neutron.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "TiaraDimensions.hh"
#include "TiaraVSourceEnergyGenerator.hh"
#include "TiaraVDirectionGenerator.hh"

TiaraPrimaryGeneratorAction::
TiaraPrimaryGeneratorAction(const TiaraVSourceEnergyGenerator& eG,
			    const TiaraVDirectionGenerator& dG,
			    const TiaraTally &tally,
			    const TiaraDimensions &tiaraDimensions) :
  fEnergyGenerator(eG.Clone()),
  fDirectionGenerator(dG.Clone()), 
  particleGun(new G4ParticleGun(1)),
  fTally(tally)
{
  particleGun->
    SetParticlePosition(G4ThreeVector(0.0, 0.0, 
				      tiaraDimensions.targetPosZ));
}

TiaraPrimaryGeneratorAction::~TiaraPrimaryGeneratorAction()
{
  delete particleGun;
  delete fDirectionGenerator;
  delete fEnergyGenerator;
}

TiaraPrimaryGeneratorAction::
TiaraPrimaryGeneratorAction(const TiaraPrimaryGeneratorAction &rhs)
  : G4VUserPrimaryGeneratorAction()
{  
  *this = rhs;
}

TiaraPrimaryGeneratorAction &
TiaraPrimaryGeneratorAction::
operator=(const TiaraPrimaryGeneratorAction &rhs)
{  
  if (this != &rhs) {
    fEnergyGenerator = rhs.fEnergyGenerator->Clone();
    fDirectionGenerator = rhs.fDirectionGenerator->Clone(); 
    particleGun = new G4ParticleGun(1);
    fTally = rhs.fTally;
  }
  return *this;
}

void TiaraPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->SetParticleDefinition(G4Neutron::NeutronDefinition());
  particleGun->SetParticleMomentumDirection(fDirectionGenerator->
					    GetDirection());  


  G4double e(fEnergyGenerator->GetEnergy());
  fTally.fill(e, 1);
  fTally.EndOfEventAction();
  particleGun->SetParticleEnergy(e);

  particleGun->GeneratePrimaryVertex(anEvent);
}

const TiaraVSourceEnergyGenerator *TiaraPrimaryGeneratorAction::
GetEnergyGenerator() const {
  return fEnergyGenerator;
}

const TiaraTally &
TiaraPrimaryGeneratorAction::GetTally() const {
  return fTally;
}
