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
// $Id: G4FinalStateChargeIncrease.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4FinalStateChargeIncrease.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateChargeIncrease::G4FinalStateChargeIncrease()
{
  lowEnergyLimit = 1 * keV;
  highEnergyLimit = 10 * MeV;

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class 4FinalStateChargeIncrease is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateChargeIncrease::~G4FinalStateChargeIncrease()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4FinalStateProduct& G4FinalStateChargeIncrease::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  product.Clear();

  //SI - Added protection against total energy deposit
  product.DoNotDepositEnergy();
  //
  product.KillPrimaryParticle();
  product.AddEnergyDeposit(0.);

  G4ParticleDefinition* definition = track.GetDefinition();
 
  G4double inK = track.GetDynamicParticle()->GetKineticEnergy();
  
  G4int finalStateIndex = cross.RandomSelect(inK,definition);

  G4int n = NumberOfFinalStates(track.GetDefinition(),finalStateIndex);
  G4double outK = inK - IncomingParticleBindingEnergyConstant(definition,finalStateIndex);

  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double electronK;
  if (definition == instance->GetIon("hydrogen")) electronK = inK*electron_mass_c2/proton_mass_c2;
  else electronK = inK*electron_mass_c2/(3728*MeV);
  
  if (outK<0)
  {
    G4String message;
    message="G4FinalStateChargeIncrease - Final kinetic energy is below 0! Process ";
  }
  
  product.AddSecondary(new G4DynamicParticle(OutgoingParticleDefinition(definition,finalStateIndex), 
					     track.GetDynamicParticle()->GetMomentumDirection(), 
					     outK));

  n = n - 1;
  
  while (n>0)
  {
    n--;
    product.AddSecondary
	(new G4DynamicParticle(G4Electron::Electron(), track.GetDynamicParticle()->GetMomentumDirection(), electronK));
  }
        
  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4FinalStateChargeIncrease::NumberOfFinalStates(G4ParticleDefinition* particleDefinition, 
						      G4int finalStateIndex )
 
{
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();
  if (particleDefinition == instance->GetIon("hydrogen")) return 2;
  if (particleDefinition == instance->GetIon("alpha+")) return 2;
  
  if (particleDefinition == instance->GetIon("helium")) 
  {    if (finalStateIndex==0) return 2;
       return 3;
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition*  G4FinalStateChargeIncrease::OutgoingParticleDefinition (G4ParticleDefinition* particleDefinition, 
									       G4int finalStateIndex) 
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  if (particleDefinition == instance->GetIon("hydrogen")) return G4Proton::Proton(); 
  if (particleDefinition == instance->GetIon("alpha+")) return instance->GetIon("alpha++");

  if (particleDefinition == instance->GetIon("alpha+")) return instance->GetIon("helium");
  {
    if (finalStateIndex==0) return instance->GetIon("alpha+"); 
    return instance->GetIon("alpha++");
  }
  return 0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4FinalStateChargeIncrease::IncomingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition, 
									   G4int finalStateIndex )
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  if(particleDefinition == instance->GetIon("hydrogen")) return 13.6*eV;
  
  if(particleDefinition == instance->GetIon("alpha+"))
  {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV
      return 54.509*eV;
  }
   
  if(particleDefinition == instance->GetIon("helium"))
  {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV

      if (finalStateIndex==0) return 24.587*eV;
      return (54.509 + 24.587)*eV;
  }  

  return 0;
}
