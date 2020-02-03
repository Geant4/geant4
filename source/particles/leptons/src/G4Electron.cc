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
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4Electron.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                         ELECTRON                               ###
// ######################################################################
G4Electron* G4Electron::theInstance = 0;

G4Electron* G4Electron::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "e-";
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

    anInstance = new G4ParticleDefinition(
                 name,  electron_mass_c2,       0.0*MeV,    -1.*eplus, 
		    1,                 0,             0,          
		    0,                 0,             0,             
	     "lepton",                 1,             0,          11,
		 true,              -1.0,          NULL,
             false,                  "e"
              );
    // Bohr Magnetron
   G4double muB =  -0.5*eplus*hbar_Planck/(electron_mass_c2/c_squared) ;
   
   anInstance->SetPDGMagneticMoment( muB * 1.00115965218076 );

  }
  theInstance = reinterpret_cast<G4Electron*>(anInstance);
  return theInstance;
}

G4Electron*  G4Electron::ElectronDefinition()
{
  return Definition();
}

G4Electron*  G4Electron::Electron()
{
  return Definition();
}


