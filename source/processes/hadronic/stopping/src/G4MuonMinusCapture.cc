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
// $Id: G4MuonMinusCapture.cc 74377 2013-10-04 08:29:54Z gcosmo $
//
//---------------------------------------------------------------------
//
// GEANT4 Class 
//
// File name:     G4MuonMinusCapture
//
// Author V.Ivanchenko 25 April 2012 
//
//
// Class Description:
//
// Base process class for stopping of mu-
//
// Modifications: 
//
//  20121003  K. Genser -- Changed the constructor argument type
//                         Used two argument base constructor
//  20121016  K. Genser -- Reverting to use one argument base c'tor
//  20121002  K. Genser -- Replaced G4MuMinusCapturePrecompound with
//                         G4CascadeInterface (Bertini)
//------------------------------------------------------------------------

#include "G4MuonMinusCapture.hh"
#include "G4HadronicProcessType.hh"
#include "G4MuonMinusBoundDecay.hh"
#include "G4HadronicInteraction.hh"
#include "G4MuonMinus.hh"
#include "G4CascadeInterface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusCapture::G4MuonMinusCapture(G4HadronicInteraction* hiptr)
  : G4HadronStoppingProcess ("muMinusCaptureAtRest")
{
  SetBoundDecay(new G4MuonMinusBoundDecay()); // Owned by InteractionRegistry
  if (!hiptr) {
    hiptr = new G4CascadeInterface(); // Owned by InteractionRegistry
  }
  RegisterMe(hiptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusCapture::~G4MuonMinusCapture()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuonMinusCapture::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4MuonMinus::MuonMinus());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusCapture::ProcessDescription(std::ostream& outFile) const
{
  outFile << "Stopping of mu- using default element selector, EM cascade"
          << " sampling and bound decay sampling.\n"
	  << "Bertini model is used for nuclear capture\n"; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
