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
// $Id: G4eScreenedRutherfordTest.cc,v 1.1 2007-05-02 17:22:39 pia Exp $
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

#include "G4eCrossSectionScreenedRutherford.hh"

int main()
{
  //  G4cout.setf( ios::scientific, ios::floatfield );

  G4eCrossSectionScreenedRutherford* cross = new G4eCrossSectionScreenedRutherford;

  // Particle definitions
  
  //  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  //  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  
  // Create a DynamicParticle  
  
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
 
  G4ParticleMomentum direction(initX,initY,initZ);
 
  G4double e0 = 0.01 * keV;
  G4double eStep = 0.0001 * keV;

  for (G4int i=0; i<1000; i++) 
    {
      G4double energy = e0 + i * eStep;
      
      G4DynamicParticle dynamicParticle(electron,direction,energy);
      
      //     dynamicParticle.DumpInfo(0);
      
      // Track 
      
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      
      G4Track track(&dynamicParticle,time,position);
      
      G4double sigma = cross->CrossSection(track);
      
      G4cout << energy/keV <<" " << sigma / cm2 << G4endl;
    }

  e0 = 0.1 * keV;
  eStep = 0.001 * keV;

  for (G4int i=0; i<901; i++) 
    {
      G4double energy = e0 + i * eStep;
      
      G4DynamicParticle dynamicParticle(electron,direction,energy);
      
      //     dynamicParticle.DumpInfo(0);
      
      // Track 
      
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      
      G4Track track(&dynamicParticle,time,position);
      
      G4double sigma = cross->CrossSection(track);
      
      G4cout << energy/keV <<" " << sigma / cm2 << G4endl;
    }

  delete cross;

  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








