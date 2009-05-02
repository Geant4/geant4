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
// $Id: printMassTable.cc,v 1.1 2009-05-02 07:23:11 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "globals.hh"
#include "tst2ParticleConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"

#include <fstream>
#include <iomanip>


int main(int ,char** ) {

  //--- open index file -----
  G4String fileName = "massTable.txt";
  std::ofstream outFile(fileName, std::ios::out );
  outFile.setf( std::ios::fixed);
  outFile.setf( std::ios::right, std::ios::adjustfield );

  
  // set readiness for the particle table
  G4ParticleTable::GetParticleTable()->SetReadiness();

  // create all particles
  tst2ParticleConstructor pConstructor;
  pConstructor.ConstructParticle();

  // loop over nuclei
  G4int iz, ia;
  G4double z,a;
  for (iz=0; iz <150; iz +=1) {
    z = iz;
    
    for (ia=iz; ia < 4*iz; ia +=1) {
      a =ia;

      // check if (A,Z) is in mass table
      if ( !G4NucleiProperties::IsInStableTable(a,z)) continue;
      
      G4double mass = G4NucleiProperties::GetNuclearMass(a, z);
      G4double massExces = G4NucleiProperties::GetMassExcess(a, z); 
      G4double bindE = G4NucleiProperties::GetBindingEnergy(a, z); 
      
      outFile  << std::setw(4) << iz << "  " 
	       << std::setw(4) << ia << "  " 
	       << std::setw(12) << massExces/keV << "   "
	       << std::setw(12) << bindE/keV << "   "
	       << std::setw(15) << mass/GeV << "   "
	       << std::endl;
    }
  }


  return EXIT_SUCCESS;
}




