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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520 
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $Id$
//
/// \file ITTrackingAction.hh 
/// \brief Implementation of the ITTrackingAction class

#include "ITTrackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4UnitsTable.hh"
//
//#include "NeuronHitCompartments.hh"
#include "G4MoleculeCounter.hh"
#include "G4MoleculeGun.hh"
#include "G4H2O.hh"
#include <G4Scheduler.hh>
#include "G4MoleculeTable.hh"
#include "math.h"
#include "Run.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITTrackingAction::ITTrackingAction(): G4UserTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITTrackingAction::~ITTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingAction::PreUserTrackingAction(const G4Track* track)
{
 // Target volumes
 G4VPhysicalVolume* volumeStep = track->GetTouchableHandle()->GetVolume();
 G4VPhysicalVolume* volumeMedium = G4PhysicalVolumeStore::GetInstance()->
                                   GetVolume("Medium");
 G4VPhysicalVolume* volumeSlice = G4PhysicalVolumeStore::GetInstance()->
                                  GetVolume("BoundingSlice"); 
 
 // count produced species in neuron
 G4Molecule* molecule = GetMolecule(track);
 //const G4String& moleculeName = molecule->GetName();
 //G4double GTime  = track->GetGlobalTime() /picosecond;
 
 Run* run = static_cast<Run*>(
            G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
 //run->MoleculeCount(moleculeName,GTime);
 
 // time steps        

 // particles outside neuron structure
 //if (volumeStep == volumeMedium || volumeStep == volumeSlice)
 //{
  //run->MoleculeCount(molecule);
 //}
 // count secondary particles in neuron
 //else //
 if (volumeStep != volumeMedium && volumeStep != volumeSlice)
 { 
  run->MoleculeCountNeuron(molecule);
      //run->MoleculeCountNeuron(moleculeName,GTime); 
     
   // number of molecules at time of 1 ns !   
   //if ( (GTime > T_1ns-3.) && (GTime <= T_1ns))
  // {  
  //  if (moleculeName =="OH^0") // OH* radical
  //  {
     //fEventAction->MoleculeCountNeuron()
  //G4MoleculeCounter::Instance()->GetNMoleculesAtTime(moleculeName, GTime);
 
 //   }
  // }
 }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
ITTrackingAction::PostUserTrackingAction(const G4Track* /*track*/)
{
}
