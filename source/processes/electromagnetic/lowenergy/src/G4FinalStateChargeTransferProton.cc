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
// $Id: G4FinalStateChargeTransferProton.cc,v 1.2 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Final state generation for charge transfer:
// incident proton hydrocarbon target
// outgoing neutral hydrogen
// Further documentation available from http://www.ge.infn.it/geant4/

// -------------------------------------------------------------------


#include "G4FinalStateChargeTransferProton.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAGenericIonsManager.hh"

G4FinalStateChargeTransferProton::G4FinalStateChargeTransferProton()
{
  name = "ChargeTransferProton";
  lowEnergyLimit = 0.1 * eV;
  highEnergyLimit = 1000. * MeV;

  // Proton binding
  protonBinding = 13.6 * eV;

  // Ionisation potentials
  G4String materialName;
  G4double potential = 0.;

  materialName = "CH";
  potential = 11.13 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "H2O";
  potential = 10.79 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CH4";
  potential = 12.51 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H6";
  potential = 11.52 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H8";
  potential = 11.08 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CH3";
  potential = 9.84  * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H5";
  potential = 8.25  * eV;
  ioniPotentialMap[materialName] = potential;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H7";
  potential = 9.10 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CH2"; 
  potential = 10.46 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H4";
  potential = 10.51 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H6";
  potential = 9.74 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CH";
  potential = 11.13  * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H3";
  potential = 9.45 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H5";
  potential = 9.90  * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H2";
  potential = 11.51 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H4";
  potential =  10.02 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C2H";
  potential =  17.42 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "Propyne";
  potential =  10.36 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "Allene";
  potential = 9.69 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H3";
  potential =  8.34 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H2";
  potential =  12.50 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "C3H";
  potential =  13.40 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "N2";
  potential =  15.58 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "O2";
  potential =  12.07 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CO2";
  potential =  13.77 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "CO";
  potential =  14.01 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "NH3";
  potential =  10.07 * eV;
  ioniPotentialMap[materialName] = potential;

  materialName = "NO";
  potential =  9.26 * eV;
  ioniPotentialMap[materialName] = potential;
}


G4FinalStateChargeTransferProton::~G4FinalStateChargeTransferProton()
{ }
 

const G4FinalStateProduct& G4FinalStateChargeTransferProton::GenerateFinalState(const G4Track& track, 
									    const G4Step& /* step */)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  G4double e = track.GetDynamicParticle()->GetKineticEnergy();

  // Get the material the track is in
  G4Material* material = track.GetMaterial();
  G4String materialName;
  materialName = material->GetName();

  // Check whether cross section data are available for the current material

  std::map<G4String,G4double,std::less<G4String> >::const_iterator pos;
  pos = ioniPotentialMap.find(materialName);

  // Generate final state only if the material is among those included in the ionisation potential
  // database; otherwise do nothing
  if (pos!= ioniPotentialMap.end())
    {
      // Ionisation potential of the medium
      G4double potential = (*pos).second;
  
      // Kinetic energy of the outgoing hydrogen ion
      G4double eProduct = e - (e * electron_mass_c2/proton_mass_c2) - potential + protonBinding;
            
      if (eProduct < 0)
	{
	  G4String message;
	  message="ChargeTransferProton::GenerateFinalState - Secondary product has negative kinetic energy";
	  G4Exception(message);
	}
      
      // Primary particle is killed
      product.KillPrimaryParticle();

      // Deposit the ionisation potential equivalent locally
      product.AddEnergyDeposit(potential);
      
      //Secondary particle
      // Create a Hydrogen ion, same direction as the incoming proton
      G4DNAGenericIonsManager* instance = G4DNAGenericIonsManager::Instance();
      G4ParticleDefinition* hydrogenDefinition = instance->GetIon("hydrogen"); 
 

      G4DynamicParticle* secondary = new G4DynamicParticle(hydrogenDefinition, 
							   track.GetDynamicParticle()->GetMomentumDirection(), 
							   eProduct);
      product.AddSecondary(secondary);
      
    }
  
  return product;
}
