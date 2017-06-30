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
/// \file exoticphysics/dmparticle/src/G4LDMPhoton.cc
/// \brief Implementation of the G4LDMPhoton class
//
// $Id: G4LDMPhoton.cc 66817 2013-01-12 16:16:08Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LDMPhoton
//  
// Description: 
//
// 15.03.17  V. Grichine based on G4Monopole
//
// Modified:
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4LDMPhoton.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

G4LDMPhoton* G4LDMPhoton::theLDMPhoton = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMPhoton::G4LDMPhoton(
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

G4LDMPhoton::~G4LDMPhoton()
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

G4LDMPhoton* G4LDMPhoton::LDMPhotonDefinition(G4double mass)
{    
  if(!theLDMPhoton) 
  {
    theLDMPhoton = new G4LDMPhoton(
       "ldmphoton",         mass,       5.317e-14*MeV,       0, 
                2,               -1,             -1,          
                0,               0,             0,             
          "boson",               0,             0,           50,
             true,            12.380*ns,             NULL);
    
    //create Decay Table
    G4DecayTable* table = new G4DecayTable();

    // create decay channels
    G4VDecayChannel** mode = new G4VDecayChannel*[3];
    // 
    mode[0] = 
      new G4PhaseSpaceDecayChannel("ldmphoton",0.999,2,"ldmhi","ldmhibar");
    // 
    mode[1] = new G4PhaseSpaceDecayChannel("ldmphoton",0.0005,2,"mu+","mu-");
    // 
    mode[2] = new G4PhaseSpaceDecayChannel("ldmphoton",0.0005,3,"pi+","pi-");

    for (G4int index=0; index <3; index++ ) table->Insert(mode[index]);

    delete [] mode;

    theLDMPhoton->SetDecayTable(table);
     
    G4cout << "LDMPhoton is created: m(GeV)= " 
           << theLDMPhoton->GetPDGMass()/GeV << G4endl;
  }
  return theLDMPhoton;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LDMPhoton* G4LDMPhoton::LDMPhoton()
{    
  return theLDMPhoton;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



