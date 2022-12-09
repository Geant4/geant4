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
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//    1st March 2007 creation by L. Desorgher based on a modification of G4Electron 	   
//    06  Nov.  2008 modified for Geant4-09-02  by Hisaya Kurashige 
//
//-------------------------------------------------------------
// 

#include "G4AdjointElectronFI.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###             ADJOINT ELECTRON FOR FORCED INTERACTION            ###
// ######################################################################
G4AdjointElectronFI* G4AdjointElectronFI::theInstance = 0;

G4AdjointElectronFI* G4AdjointElectronFI::Definition()
{ 
  
  if (theInstance !=0) return theInstance;
  const G4String name = "adj_e-_FI";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding

  // use constants in CLHEP
  //  static const double electron_mass_c2 = 0.51099906 * MeV;

  // NOTE  : electric charge and magnetic moment is opposite 
  //         compared with real electron, because adjoint particles
  //         are used to back-trace

   anInstance = new G4ParticleDefinition(
                 name,  electron_mass_c2,       0.0*MeV,     1.*eplus, 
		    1,                 0,             0,          
		    0,                 0,             0,             
	     "adjoint",                1,             0,     90000011,
	 	  true,             -1.0,          NULL,
                 false,           "adj_lepton"
              );
   // Bohr Magnetron
   G4double muB =  0.5*eplus*hbar_Planck/(electron_mass_c2/c_squared) ;
   anInstance->SetPDGMagneticMoment( muB * 2.* 1.0011596521859 );

  }
  theInstance = static_cast<G4AdjointElectronFI*>(anInstance);
  return theInstance;
}

G4AdjointElectronFI*  G4AdjointElectronFI::AdjointElectronFIDefinition()
{
  return Definition();
}

G4AdjointElectronFI*  G4AdjointElectronFI::AdjointElectronFI()
{
  return Definition();
}


