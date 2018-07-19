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
// $Id: G4EmCaptureCascade.cc 101422 2016-11-17 10:41:23Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class file 
//
// File name:  G4EmCaptureCascade
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 22 April 2012 on base of G4MuMinusCaptureCascade
//
//
//-----------------------------------------------------------------------------
//
// Modifications: 
//
//-----------------------------------------------------------------------------

#include "G4EmCaptureCascade.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh" 
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4NucleiProperties.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCaptureCascade::G4EmCaptureCascade()
  : G4HadronicInteraction("emCaptureCascade")
{ 
  theElectron = G4Electron::Electron();
  theGamma = G4Gamma::Gamma();
  fMuMass = G4MuonMinus::MuonMinus()->GetPDGMass();
  fTime = 0.0;

  // Calculate the Energy of K Mesoatom Level for this Element using
  // the Energy of Hydrogen Atom taken into account finite size of the
  // nucleus 
  static const G4int nlevels = 28;
  static const G4int listK[nlevels] = {
      1, 2,  4,  6,  8, 11, 14, 17, 18, 21, 24,
     26, 29, 32, 38, 40, 41, 44, 49, 53, 55,
     60, 65, 70, 75, 81, 85, 92};
  static const G4double listKEnergy[nlevels] = {
     0.00275, 0.011, 0.043, 0.098, 0.173, 0.326,
     0.524, 0.765, 0.853, 1.146, 1.472,
     1.708, 2.081, 2.475, 3.323, 3.627, 
     3.779, 4.237, 5.016, 5.647, 5.966,
     6.793, 7.602, 8.421, 9.249, 10.222,
    10.923,11.984};

  fKLevelEnergy[0] = 0.0;
  fKLevelEnergy[1] = listKEnergy[0];
  G4int idx = 1;
  for(G4int i=1; i<nlevels; ++i) {
    G4int z1 = listK[idx];
    G4int z2 = listK[i];
    if(z1+1 < z2) {
      G4double dz = G4double(z2 - z1);
      G4double y1 = listKEnergy[idx]/G4double(z1*z1);
      G4double y2 = listKEnergy[i]/G4double(z2*z2);
      for(G4int z=z1+1; z<z2; ++z) {
        fKLevelEnergy[z] = (y1 + (y2 - y1)*(z - z1)/dz)*z*z;
      }
    }
    fKLevelEnergy[z2] = listKEnergy[i];
    idx = i;  
  }
  for(G4int i = 0; i<14; ++i) { fLevelEnergy[i] = 0.0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCaptureCascade::~G4EmCaptureCascade()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadFinalState* 
G4EmCaptureCascade::ApplyYourself(const G4HadProjectile& projectile, 
				  G4Nucleus& targetNucleus)
{
  result.Clear();
  result.SetStatusChange(isAlive);
  fTime = projectile.GetGlobalTime();

  G4int Z = targetNucleus.GetZ_asInt(); 
  G4int A = targetNucleus.GetA_asInt(); 
  G4double massA = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double mass = fMuMass * massA / (fMuMass + massA) ;
  G4double e = 13.6 * eV * (Z * Z) * mass/ electron_mass_c2;

  // precise corrections of energy only for K-shell
  fLevelEnergy[0] = fKLevelEnergy[std::min(Z, 92)];
  for(G4int i=1; i<14; ++i) {
    fLevelEnergy[i] = e/(G4double)((i+1)*(i+1));
  }

  G4int nElec  = Z;
  G4int nAuger = 1;
  G4int nLevel = 13;
  G4double pGamma = (Z*Z*Z*Z);

  // Capture on 14-th level
  G4double edep = fLevelEnergy[13];
  AddNewParticle(theElectron,edep);
  G4double deltaE;

  // Emit new photon or electron
  // Simplified model for probabilities
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  do {

    // case of Auger electrons
    if((nAuger < nElec) && ((pGamma + 10000.0) * G4UniformRand() < 10000.0) ) {
      ++nAuger;
      deltaE =  fLevelEnergy[nLevel-1] -  fLevelEnergy[nLevel];
      --nLevel;
      AddNewParticle(theElectron, deltaE);

    } else {

      // Case of photon cascade, probabilities from
      // C.S.Wu and L.Wilets, Ann. Rev. Nuclear Sci. 19 (1969) 527.

      G4double var = (10.0 + G4double(nLevel - 1) ) * G4UniformRand();
      G4int iLevel = nLevel - 1 ;
      if(var > 10.0) iLevel -= G4int(var-10.0) + 1;
      if( iLevel < 0 ) iLevel = 0;
      deltaE =  fLevelEnergy[iLevel] -  fLevelEnergy[nLevel];
      nLevel = iLevel;
      AddNewParticle(theGamma, deltaE);
    }
    edep += deltaE;

    // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
  } while( nLevel > 0 );

  result.SetLocalEnergyDeposit(edep);
  return &result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCaptureCascade::ModelDescription(std::ostream& outFile) const
{
  outFile << "Simulation of electromagnetic cascade from capture level"
	  << " to K-shell of the mesonic atom\n."
	  << "Probabilities of gamma and Auger transitions from\n"
	  << "  N.C.Mukhopadhyay Phys. Rep. 30 (1977) 1.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
