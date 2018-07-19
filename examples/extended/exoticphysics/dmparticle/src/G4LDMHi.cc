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
/// \file exoticphysics/dmparticle/src/G4LDMHi.cc
/// \brief Implementation of the G4LDMHi class
//
// $Id: G4LDMHi.cc 66817 2013-01-12 16:16:08Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LDMHi
//  
// Description: 
//
// 15.03.17  V. Grichine based on G4Monopole
//
// Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4LDMHi.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMHi* G4LDMHi::theLDMHi = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMHi::G4LDMHi(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable)
 : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity,
           iConjugation, iIsospin, iIsospin3, gParity, pType,
           lepton, baryon, encoding, stable, lifetime, decaytable )
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMHi::~G4LDMHi()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//     
//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table 
//
//

G4LDMHi* G4LDMHi::LDMHiDefinition(G4double mass)
{    
  if(!theLDMHi) 
  {
    theLDMHi = new G4LDMHi(
       "ldmhi",         mass,       0.0*MeV,       0, 
                0,               -1,             -1,          
                0,               0,             0,             
          "boson",               0,             0,           51,
             true,            -1.0,             0);
    
    
    G4cout << "LDMHi is created: m(GeV)= " << theLDMHi->GetPDGMass()/GeV 
           << G4endl;
  }
  return theLDMHi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMHi* G4LDMHi::LDMHi()
{    
  return theLDMHi;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


