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
// $Id: G4FinalStateChargeDecrease.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4FinalStateChargeDecrease.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateChargeDecrease::G4FinalStateChargeDecrease()
{
  lowEnergyLimit = 1 * keV;
  highEnergyLimit = 10 * MeV;

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4FinalStateChargeDecrease is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateChargeDecrease::~G4FinalStateChargeDecrease()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4FinalStateProduct& G4FinalStateChargeDecrease::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  product.Clear();

  G4double inK = track.GetDynamicParticle()->GetKineticEnergy();

  G4ParticleDefinition* definition = track.GetDefinition();
  G4int finalStateIndex = cross.RandomSelect(inK,definition);
 
  G4int n = NumberOfFinalStates(definition, finalStateIndex);
  G4double waterBindingEnergy = WaterBindingEnergyConstant(definition, finalStateIndex);
  G4double outgoingParticleBindingEnergy = OutgoingParticleBindingEnergyConstant(definition, finalStateIndex);
  
  G4double outK = 0.;
  if (definition==G4Proton::Proton())
    outK = inK - n*(inK*electron_mass_c2/proton_mass_c2) - waterBindingEnergy + outgoingParticleBindingEnergy;
  else
    outK = inK - n*(inK*electron_mass_c2/(3728*MeV)) - waterBindingEnergy + outgoingParticleBindingEnergy;
  
  if (outK<0)
  {
    G4String message;
    message="ChargeDecreaseDingfelder::GenerateFinalState - Final kinetic energy is below 0! Process ";
    G4Exception(message);
  }
  
  //SI - Added protection against total energy deposit
  product.DoNotDepositEnergy();
  //
  product.KillPrimaryParticle();

  product.AddEnergyDeposit(waterBindingEnergy);

  G4DynamicParticle* aSecondary = new G4DynamicParticle(OutgoingParticleDefinition(definition, finalStateIndex), 
							track.GetDynamicParticle()->GetMomentumDirection(), 
							outK);
  
  product.AddSecondary(aSecondary);
  
  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4FinalStateChargeDecrease::NumberOfFinalStates(G4ParticleDefinition* particleDefinition, 
						      G4int finalStateIndex )
  
{
  if (particleDefinition == G4Proton::Proton()) return 1;
  
  G4DNAGenericIonsManager*instance;
  instance = G4DNAGenericIonsManager::Instance();
  
  if (particleDefinition == instance->GetIon("alpha++") ) 
  {
    if (finalStateIndex==0)  return 1;
    return 2;   
  }
  
  if (particleDefinition == instance->GetIon("alpha+") ) return 1;
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4FinalStateChargeDecrease::OutgoingParticleDefinition (G4ParticleDefinition* particleDefinition, 
									      G4int finalStateIndex) 
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  
  if (particleDefinition == G4Proton::Proton()) return instance->GetIon("hydrogen"); 
  
  if (particleDefinition == instance->GetIon("alpha++") ) 
  {
    if (finalStateIndex == 0) return instance->GetIon("alpha+"); 
    return instance->GetIon("helium");    
  }
  
  if (particleDefinition == instance->GetIon("alpha+") ) return instance->GetIon("helium");    
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateChargeDecrease::WaterBindingEnergyConstant(G4ParticleDefinition* particleDefinition, 
								G4int finalStateIndex )
{
  // Ionization energy of first water shell
  // Rad. Phys. Chem. 59 p.267 by Dingf. et al.
  // W + 10.79 eV -> W+ + e-
  
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  
  if(particleDefinition == G4Proton::Proton()) return 10.79*eV;

  if (particleDefinition == instance->GetIon("alpha++") ) 
  {
      // Binding energy for    W+ -> W++ + e-    10.79 eV
      // Binding energy for    W  -> W+  + e-    10.79 eV
      //
      // Ionization energy of first water shell
      // Rad. Phys. Chem. 59 p.267 by Dingf. et al.

      if (finalStateIndex == 0) return 10.79*eV;
  
      return 10.79*2*eV;
  }

  if (particleDefinition == instance->GetIon("alpha+") ) 
  {
      // Binding energy for    W+ -> W++ + e-    10.79 eV
      // Binding energy for    W  -> W+  + e-    10.79 eV
      //
      // Ionization energy of first water shell
      // Rad. Phys. Chem. 59 p.267 by Dingf. et al.

      return 10.79*eV;
  }
  
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateChargeDecrease::OutgoingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition, 
									   G4int finalStateIndex)
{
  G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
  
  if(particleDefinition == G4Proton::Proton()) return 13.6*eV;

  if (particleDefinition == instance->GetIon("alpha++") ) 
  {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV

      if (finalStateIndex==0)	return 54.509*eV;
      
      return (54.509 + 24.587)*eV;
  }
  
  if (particleDefinition == instance->GetIon("alpha+") ) 
  {
      // Binding energy for    He+ -> He++ + e-    54.509 eV
      // Binding energy for    He  -> He+  + e-    24.587 eV
      
      return 24.587*eV;
  }  
  
  return 0;
}

