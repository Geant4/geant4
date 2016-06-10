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
// $Id: G4StopElementSelector.cc 88066 2015-01-27 10:41:12Z gcosmo $
//
// File: G4StopElementSelector
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 2 April 2000
//
// Modifications: 
// 18/08/2000  V.Ivanchenko Update description
// 17/05/2006  V.Ivanchenko Cleanup
// 02/10/2007  V.Ivanchenko Fixed typo in computation of Lambda-factor
//                          proposed by Victor Pec 
// 04/23/2013  K.Genser     used new G4MuonMinusBoundDecay
//                          in GetMuonCaptureRate and GetMuonDecayRate
//---------------------------------------------------------------------

#include "G4StopElementSelector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" 
#include "G4Material.hh"
#include "G4MuonMinusBoundDecay.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// constructor
G4StopElementSelector::G4StopElementSelector()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// destructor
G4StopElementSelector::~G4StopElementSelector()
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Element* G4StopElementSelector::GetElement(const G4Material* aMaterial)
{
  // Fermi-Teller Z-low of mu- capture and exceptions 
  // for halogens and oxigen.
  // N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.
  G4int i;
  G4double Z;
  const G4int numberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();

  if(1 == numberOfElements) return (*theElementVector)[0];
    
  const G4double* theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();

  G4double sum = 0.0;
  for ( i=0; i < numberOfElements; i++ ) {

    Z = (*theElementVector)[i]->GetZ();

      // Halogens
    if( (9.0 == Z) || (17.0 == Z) || (35.0 == Z) || (53.0 == Z) || (85.0 == Z) ) {
      sum += 0.66 * Z * theAtomicNumberDensity[i] ; 

      // Oxigen
    } else if( 8.0 == Z ) {
      sum += 0.56 * Z * theAtomicNumberDensity[i] ; 

      // Others
    } else {
      sum +=        Z * theAtomicNumberDensity[i] ; 
    }
  }

  G4double random = G4UniformRand() * sum;
  sum = 0.0 ;
  i   = -1;

  // Selection of element
  do {
    i++;
    Z = (*theElementVector)[i]->GetZ();

      // Galogens
    if( (9.0 == Z) || (17.0 == Z) || (35.0 == Z) || (53.0 == Z) || (85.0 == Z) ) {
      sum += 0.66 * Z * theAtomicNumberDensity[i] ; 

      // Oxigen
    } else if( 8.0 == Z ) {
      sum += 0.56 * Z * theAtomicNumberDensity[i] ; 

      // Others
    } else {
      sum +=        Z * theAtomicNumberDensity[i] ; 
    }
  } while ( (sum < random) && (i < numberOfElements - 1) );

  return (*theElementVector)[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double  G4StopElementSelector::GetMuonCaptureRate(G4double Z, G4double A)
{
  return G4MuonMinusBoundDecay::GetMuonCaptureRate(G4int(Z),G4int(A));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double  G4StopElementSelector::GetMuonDecayRate(G4double Z, G4double /* A */)
{
  return G4MuonMinusBoundDecay::GetMuonDecayRate(G4int(Z));
}
