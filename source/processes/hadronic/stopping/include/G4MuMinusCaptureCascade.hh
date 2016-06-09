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
// $Id: G4MuMinusCaptureCascade.hh,v 1.10 2007-07-05 18:19:14 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   G4MuMinusCaptureCascade physics process --------
//   Vladimir Ivanchenko, April 2000
//
// Modified:
// 14.11.06 Add inline functions (V.Ivanchenko)
//
//-----------------------------------------------------------------------------

#ifndef G4MuMinusCaptureCascade_h
#define G4MuMinusCaptureCascade_h 1

#include "globals.hh"
#include "Randomize.hh" 
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"

class G4GHEKinematicsVector;

class G4MuMinusCaptureCascade 
{ 
public:
 
  G4MuMinusCaptureCascade();
 
  ~G4MuMinusCaptureCascade();

  G4int DoCascade(const G4double Z, const G4double A,
                  G4GHEKinematicsVector* Cascade);

  void DoBoundMuonMinusDecay(G4double Z, 
			     G4int* nCascade, G4GHEKinematicsVector* Cascade);

  G4double GetKShellEnergy(G4double Z);

  G4ThreeVector& GetRandomVec();

private:

  G4double GetLinApprox(G4int N, const G4double* X, const G4double* Y, 
			G4double Xuser);

  void AddNewParticle(G4ParticleDefinition* aParticle,
		      G4ThreeVector& Momentum,
		      G4double mass,
		      G4int* nParticle,
		      G4GHEKinematicsVector* Cascade);

  // hide assignment operator as private 
  G4MuMinusCaptureCascade& operator=(const G4MuMinusCaptureCascade &right);
  G4MuMinusCaptureCascade(const G4MuMinusCaptureCascade& );

  G4double Emass, MuMass;
  G4ParticleDefinition* theElectron;
  G4ParticleDefinition* theGamma;
  G4ThreeVector randomVect;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuMinusCaptureCascade::GetLinApprox(G4int N, 
						      const G4double* X, 
						      const G4double* Y, 
						      G4double Xuser)
{
  G4double Yuser;
  if(Xuser <= X[0])        Yuser = Y[0];
  else if(Xuser >= X[N-1]) Yuser = Y[N-1];
  else {
    G4int i;
    for (i=1; i<N; i++){
      if(Xuser <= X[i]) break; 
    }    

    if(Xuser == X[i]) Yuser = Y[i];
    else Yuser = Y[i-1] + (Y[i] - Y[i-1])*(Xuser - X[i-1])/(X[i] - X[i-1]);
  }
  return Yuser;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4ThreeVector& G4MuMinusCaptureCascade::GetRandomVec() 
{
  //
  // generate uniform vector
  //
  G4double cost = 2.0 * G4UniformRand() - 1.0;
  G4double sint = std::sqrt((1.0 - cost)*(1.0 + cost));
  G4double Phi  = twopi * G4UniformRand();
  randomVect = G4ThreeVector(sint * std::cos(Phi), sint * std::sin(Phi), cost);
  return randomVect;
}

#endif


