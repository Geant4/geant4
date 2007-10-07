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
// $Id: G4FinalStateTest.cc,v 1.1 2007-10-07 12:59:37 pia Exp $
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
//#include "G4ParticleTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

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
    
    // Track 
    G4ThreeVector position(0.,0.,0.);
    G4double time = 0. ;
    const G4Track track(&dynamicParticle,time,position);
    
    // Step (dummy)
    const G4Step step;

    // Generate final state
    const G4FinalStateProduct product = finalState->GenerateFinalState(track,step);
    
    G4double energyDeposit = product.GetEnergyDeposit();
    G4int nSecondaries = product.NumberOfSecondaries();
      
    G4cout << energyDeposit/keV <<" keV deposited" << G4endl;
    G4cout << nSecondaries <<" secondaries generated" << G4endl;
 
   
   delete finalState;
  
   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








