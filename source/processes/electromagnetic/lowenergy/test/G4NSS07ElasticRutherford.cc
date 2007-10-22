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
// $Id: G4NSS07ElasticRutherford.cc,v 1.1 2007-10-22 10:21:25 pia Exp $
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

#include "G4CrossSectionElasticScreenedRutherford.hh"
#include "G4CrossSectionExcitationEmfietzoglou.hh"
#include "G4CrossSectionExcitationEmfietzoglouPartial.hh"
#include "G4CrossSectionExcitationBorn.hh"
#include "G4CrossSectionExcitationBornPartial.hh"
#include "G4CrossSectionIonisationBornElectron.hh"
#include "G4CrossSectionIonisationBorn.hh"

int main()
{
  //  G4cout.setf( ios::scientific, ios::floatfield );

  G4CrossSectionElasticScreenedRutherford* cross = new G4CrossSectionElasticScreenedRutherford;
  // G4CrossSectionExcitationEmfietzoglou* cross = new G4CrossSectionExcitationEmfietzoglou;
  // G4CrossSectionExcitationEmfietzoglouPartial* cross = new G4CrossSectionExcitationEmfietzoglouPartial;
  // G4CrossSectionExcitationBornPartial* cross = new G4CrossSectionExcitationBornPartial;
  // G4CrossSectionIonisationBornElectron* cross = new G4CrossSectionIonisationBornElectron;

  // G4CrossSectionIonisationBorn* cross = new G4CrossSectionIonisationBorn;

  // Particle definitions
  
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  //  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  //   G4ParticleDefinition* positron = G4Positron::PositronDefinition();
 
  // Create a DynamicParticle  
  
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
 
  G4ParticleMomentum direction(initX,initY,initZ);

  std::vector<G4double> mason;

  // Itikawa & Mason, J. Phys. Chem. Ref. Data 34 (1), pp. 1-22, 2005, table 4. 
  mason.push_back(1);
  mason.push_back(2);
  mason.push_back(4);
  mason.push_back(6);
  mason.push_back(10);
  mason.push_back(20);
  mason.push_back(30);
  mason.push_back(40);
  mason.push_back(50);
  mason.push_back(60);
  mason.push_back(70);
  mason.push_back(80);
  mason.push_back(90);
  mason.push_back(100);

 G4int nE = mason.size();


  for (G4int i=0; i<nE; i++) 
    {
      G4double energy = mason[i] * eV;
      
      G4DynamicParticle dynamicParticle(electron,direction,energy);
      //      G4DynamicParticle dynamicParticle(proton,direction,energy);
      //      G4DynamicParticle dynamicParticle(positron,direction,energy);
       
      //     dynamicParticle.DumpInfo(0);
     
      // Track 
      
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      
      G4Track track(&dynamicParticle,time,position);
      
      G4double sigma = cross->CrossSection(track);
      // G4double sigma = cross->CrossSection(energy,0);
      // G4double sigma = cross->CrossSection(energy,electron::GetParticleName());

      // G4cout << energy/eV <<" eV, cross section = " << sigma / (cm*cm) << " cm-2" << G4endl;

      G4cout << energy/eV <<" " << sigma / cm2 << G4endl;

      // G4int level = cross->RandomSelect(energy);

      //   G4cout << "Level = " << level << G4endl;
    }

  delete cross;

  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}








