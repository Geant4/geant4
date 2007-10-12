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
// $Id: G4FinalStateTest.cc,v 1.2 2007-10-12 16:39:46 pia Exp $
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


#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4FinalStateProduct.hh"
#include "G4DummyFinalState.hh"


int main()
{
  //  G4cout.setf( ios::scientific, ios::floatfield );

  G4DummyFinalState* finalState = new G4DummyFinalState;

  // Create particle track

  // Particle definition
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  
  // Create a DynamicParticle  
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  G4ParticleMomentum direction(initX,initY,initZ);
  
  G4cout << "Enter energy " << G4endl;
  G4double energy;
  G4cin >> energy;

  energy = energy * keV;
 
  G4DynamicParticle dynamicParticle(electron,direction,energy);
    
    //     dynamicParticle.DumpInfo(0);    

  // Materials
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  // Dump the material table
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*materialTable)[mat]->GetName()
	     << G4endl;
    }

  G4double dimX = 1 * mm;
  G4double dimY = 1 * mm;
  G4double dimZ = 1 * mm;
  
  // Geometry 

  G4Box* theFrame = new G4Box ("Frame",dimX, dimY, dimZ);
  G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame,
						      water,
						      "LFrame", 0, 0, 0);
  logicalFrame->SetMaterial(water); 
  G4PVPlacement* physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",logicalFrame,0,false,0);
  

  // Track 
  G4ThreeVector position(0.,0.,0.);
  G4double time = 0. ;
  
  G4Track* track = new G4Track(&dynamicParticle,time,position);
  // Do I really need this?
  G4GRSVolume* touche = new G4GRSVolume(physicalFrame, 0, position);   
  //  track->SetTouchable(touche);
  
  G4Step* step = new G4Step();  
  step->SetTrack(track);
  
  G4StepPoint* point = new G4StepPoint();
  point->SetPosition(position);
  point->SetMaterial(water);
  G4double safety = 10000.*cm;
  point->SetSafety(safety);
  step->SetPreStepPoint(point);
  
  G4StepPoint* newPoint = new G4StepPoint();
  G4ThreeVector newPosition(0.,0.,0.05*mm);
  newPoint->SetPosition(newPosition);
  newPoint->SetMaterial(water);
  step->SetPostStepPoint(newPoint);
  
  step->SetStepLength(1*micrometer);
  
  track->SetStep(step); 
  
  
  // Generate final state
  const G4FinalStateProduct product = finalState->GenerateFinalState(*track,*step);
  
  G4double energyDeposit = product.GetEnergyDeposit();
  G4int nSecondaries = product.NumberOfSecondaries();
  
  G4cout << energyDeposit/keV <<" keV deposited" << G4endl;
  G4cout << nSecondaries <<" secondaries generated" << G4endl;

  std::vector<G4DynamicParticle*> secondaries = product.GetSecondaries();
  
  for (G4int i = 0; i < nSecondaries; i++) 
    {  
      G4DynamicParticle* finalParticle = secondaries[i];
      
      G4double e  = finalParticle->GetTotalEnergy();
      G4double eKin = finalParticle->GetKineticEnergy();
      G4double px = (finalParticle->GetMomentum()).x();
      G4double py = (finalParticle->GetMomentum()).y();
      G4double pz = (finalParticle->GetMomentum()).z();
      G4String particleName = finalParticle->GetDefinition()->GetParticleName();
      G4cout  << "==== Final " 
	      <<  particleName  << " "  
	      << "energy: " <<  e/keV  << " keV,  " 
	      << "eKin: " <<  eKin/keV  << " keV, " 
	      << "(px,py,pz): ("
	      <<  px/keV  << "," 
	      <<  py/keV  << ","
	      <<  pz/keV  << ") keV "
	      <<  G4endl;     
    }
  
  
  
  delete finalState;
  
  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








