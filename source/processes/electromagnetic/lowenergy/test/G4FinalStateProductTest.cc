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
// $Id: G4FinalStateProductTest.cc,v 1.2 2007-10-12 16:39:46 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
///
// -------------------------------------------------------------------
//      Author:        Maria Grazia Pia
// 
//      Creation date: 6 August 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <vector>


#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
//#include "G4ParticleTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4FinalStateProduct.hh"

int main()
{
  //  G4cout.setf( ios::scientific, ios::floatfield );

  G4FinalStateProduct* product = new G4FinalStateProduct;

  // Dump initial information

  G4double deposit = product->GetEnergyDeposit();
  G4int nSecondaries = product->NumberOfSecondaries();

  G4cout << "Initial product: deposit = " 
	 << deposit
	 << ", "
	 << nSecondaries
	 << " secondaries"
	 <<G4endl;
  G4bool kill = product->PrimaryParticleIsKilled();

  if (kill)
    {
      G4cout << "Kill incident particle" << G4endl;
    }
  else
    {
      G4cout << "Do not kill incident particle" << G4endl;    
    }

  // Particle definitions
  //  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  
  // Create a DynamicParticle  

  G4cout << "Enter energy " << G4endl;
  G4double energy;
  G4cin >> energy;

  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  G4ParticleMomentum direction(initX,initY,initZ);
  
  G4DynamicParticle dynamicParticle(electron,direction,energy);
  
  //     dynamicParticle.DumpInfo(0);
  
  // Reset empty final state (make sure Clear is harmless)
  product->Clear();
  
  // Add secondary particle and local energy deposit (10% of particle energy)

  G4double localEnergy = 0.1 * energy;
  product->AddEnergyDeposit(localEnergy);
  product->AddSecondary(&dynamicParticle);
  product->KillIncidentParticle();

  // Retrieve final state products

  G4int nProducts = product->NumberOfSecondaries();
  G4double productDeposit = product->GetEnergyDeposit();
  std::vector<G4DynamicParticle*> products = product->GetSecondaries();

  G4DynamicParticle* product0 = products[0]; 
  if (!products.empty())
    {
      //   product0 = products[0];
    }

  if (product0 != 0)
    {
      G4double charge = product0->GetCharge();
      G4double eKin = product0->GetKineticEnergy();
      G4double mass = product0->GetMass();
      
      G4cout << nProducts << " secondary products - "
	     << "Charge = " << charge 
	     << ", kinetic E = " << eKin 
	     << ", mass = " << mass 
	     << G4endl
	     << "Local energy deposit = " << productDeposit
	     << G4endl; 
    }

  kill = product->PrimaryParticleIsKilled();

  if (kill)
    {
      G4cout << "Kill incident particle after production" << G4endl;
    }
  else
    {
      G4cout << "Do not kill incident particle after production" << G4endl;    
    }

 // Clear final state produced
  product->Clear();

  nProducts = product->NumberOfSecondaries();
  productDeposit = product->GetEnergyDeposit();
  G4cout << nProducts << " secondary products  after Clear - "
	 << "Local energy deposit after Clear = " << productDeposit
	 << G4endl; 
  
 std::vector<G4DynamicParticle*> products2 = product->GetSecondaries();

 if (products2.empty())
   {
     G4cout << "Vector is empty after Clear" << G4endl;
   }

  if (product0 != 0)
    {
      G4double charge = product0->GetCharge();
      G4double eKin = product0->GetKineticEnergy();
      G4double mass = product0->GetMass();
      
      G4cout << "Old product is still alive - "
	     << "Charge = " << charge 
	     << ", kinetic E = " << eKin 
	     << ", mass = " << mass 
	     << G4endl; 
      delete product0;

    }

  delete product;

  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








